import pprint
import subprocess
import config

from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

try:
    from gt.core import *
    from gt.extended import *
    from gt.annotationsketch import *
    from gt.annotationsketch.custom_track import CustomTrack
    from gt.core.gtrange import Range
except ImportError:
    print 'Warning: Could not import gt, sketch needs to be called externally!'

#a contains b as subinterval
def contains(a,b):
    return a[0] <= b[0] and a[1] >= b[1]

#remove a from b
def remove(a,b):
    return [[b[0],a[0]], [a[1], b[1]]]

#get intervals of a cropped by b
def crop(a, b):
    if not len(b):
        return [a]

    #sort children
    b.sort(key=lambda x: x[0])

    cropped = [[a[0], b[0][0]]]
    for i in range(len(b) - 1):
        cropped.append([b[i][1], b[i+1][0]])
    cropped.append([b[-1][1], a[1]])

    return cropped

def sketch(nested, genes):
    print 'Running GT annotation sketch...'
    for g in genes:
        createGFF(genes[g], nested[g])
        visualize(g)

def createGFF(gene, nested):
    #find closest parent
    if nested:
        nested[-1]['parent'] = -1
    for i in range(len(nested) - 1):
        nested[i]['parent'] = -1
        for j in range(i+1, len(nested)):
            if contains(nested[j]['location'], nested[i]['location']):
                nested[i]['parent'] = j
                break

    #append direct children  
    for i in reversed(range(len(nested))):
        nested[i]['direct_children_locations'] = []        
        parent = nested[i]['parent']
        if parent != -1:
            nested[parent]['direct_children_locations'].append(nested[i]['location'])
          
    #GFF
    rec = SeqRecord(gene['sequence'], gene['id'])
    features = []

    for i in range(len(nested)):
        #insert BaseLine
        features.append(SeqFeature(FeatureLocation(nested[i]['location'][0], nested[i]['location'][1]),
                    type='te_base', 
                    strand=0, 
                    qualifiers={'ID': 'TE_BASE {}'.format(i)}))

        #Insert croppedIntervals
        children = nested[i]['direct_children_locations']
        cropped = crop(nested[i]['location'], children)
        for subinterval in cropped:
            features.append(SeqFeature(FeatureLocation(subinterval[0],subinterval[1]), 
                                        type='te', 
                                        strand=0,
                                        qualifiers={'ID': 'TE {}'.format(i), 'Parent': 'TE_BASE {}'.format(i)}))

        #Insert elements (domains, ppt, pbs, ...)
        j = 0
        for domain in nested[i]['features']['domains']:
            sign = (lambda x: x and (1, -1)[x < 0])(domain['frame'][0])
            if sign < 0:
                domain['location'] = [domain['location'][1], domain['location'][0]]
            
            overlap = [x for x in children if contains(domain['location'], x)]

            cropped_domain = crop(domain['location'], overlap)
            for part in cropped_domain:
                features.append(SeqFeature(FeatureLocation(part[0], part[1]),
                                                    type=domain['type'],
                                                    strand=sign,
                                                    #qualifiers={'Parent': 'DOMAIN {}-{}'.format(i,j)}))
                                                    qualifiers={'ID': 'DOMAIN {}-{}'.format(i,j),'Parent': 'TE_BASE {}'.format(i)}))

            j += 1
        
    rec.features = features

    #possible other formats as well
    with open('data/' + gene['id'] + '.gff', 'w+') as out_handle:
        gff = GFF.write([rec], out_handle)
        

def visualize(gene_id):
    if config.call_gt_sketch_externally:
        process = subprocess.Popen([config.gt_sketch_path] + config.gt_sketch_args + ['data/{}.png'.format(gene_id)] + ['data/{}.gff'.format(gene_id)])    
        return

    #Visualize using python
    #Style
    style = Style()
    style.load_file('gt.style')

    #Feature index
    feature_index = FeatureIndexMemory()

    #add gff file
    feature_index.add_gff3file('data/' + gene_id + '.gff')

    #create diagram for first sequence ID
    seqid = feature_index.get_first_seqid()
    range = feature_index.get_range_for_seqid(seqid)
    diagram = Diagram.from_index(feature_index, seqid, range, style)

    #create layout
    layout = Layout(diagram, 1400, style)
    height = layout.get_height()

    #create canvas
    canvas = CanvasCairoFile(style, 1400, height)
    #canvas = CanvasCairoFilePDF(style, 1400, height)

    #sketch layout on canvas
    layout.sketch(canvas)

    #write to file
    pngfile = 'data/' + gene_id + '.png'
    canvas.to_file(pngfile)
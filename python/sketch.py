import math
import pprint

from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from gt.core import *
from gt.extended import *
from gt.annotationsketch import *
from gt.annotationsketch.custom_track import CustomTrack
from gt.core.gtrange import Range

def containsSubinterval(a, b):
    #a contains b as subinterval
    return b[0] >= a[0] and b[1] <= a[1] 

def removeSubinterval(a,b):
    #remove subinterval b from a
    return [[a[0], b[0]], [b[1], a[1]]]

def treeToGFF(gene, tree):
    #crop intervals
    croppedTree = []
    for i in range(len(tree)):
        #find direct children
        directChildren = []
        for j in reversed(range(0,i)):
            directChild = True
            for dc in directChildren:
                if containsSubinterval(dc, tree[j]):
                    directChild = False
            if not directChild:
                continue
            if containsSubinterval(tree[i], tree[j]):
                directChildren.append(tree[j])
        directChildren.sort()
        node = [tree[i]]
        for child in directChildren:
            node = node[:-1] + removeSubinterval(node[-1], child)
        croppedTree.append(node)
    
    #GFF
    rec = SeqRecord(gene['sequence'], gene['id'])
    features = []
    i = 0
    for interval in croppedTree:
        #insert base
        feat = SeqFeature(FeatureLocation(interval[0][0], interval[-1][1]), type='gene', strand=0, qualifiers={'ID': 'Transposon {}'.format(i)})
        feat.sub_features = []
        #insert all intervals
        for subinterval in interval:
            feat.sub_features.append(SeqFeature(FeatureLocation(subinterval[0],subinterval[1]), type='te_base', strand=0))

        features.append(feat)
        i += 1

    rec.features = features
    #possible other formats as well
    with open('data/' + gene['id'] + '.gff', 'w+') as out_handle:
        gff = GFF.write([rec], out_handle)

def sketch(genes):
    print 'Running GT annotation sketch...'
    for g in genes:
        treeToGFF(genes[g], genes[g]['nested'])
        visualize(genes[g])

def geneToGFF(gene):
    #insert gene
    rec = SeqRecord(gene['sequence'], gene['id'])
    qualifiers = {'ID': 'Gene'}
    features = []
    features.append(SeqFeature(FeatureLocation(0, len(gene['sequence'])), type='gene', strand=0, qualifiers=qualifiers))

    #insert LTR finder results
    for transposon in gene['ltr_finder']:
        LTR_feature = SeqFeature(FeatureLocation(transposon['location'][0], transposon['location'][1]), 
            type='LTR_finder_TE', strand=1)
        LTR_feature.sub_features = []
        if not math.isnan(transposon['ppt'][0]):
            strand = 1
            if transposon['ppt'][0] > transposon['ppt'][1]:
                strand = -1
            LTR_feature.sub_features.append(SeqFeature(FeatureLocation(transposon['ppt'][0], transposon['ppt'][1]), 
                type='LTR_ppt', strand=1))
        if not math.isnan(transposon['pbs'][0]):
            strand = 1
            if transposon['pbs'][0] > transposon['pbs'][1]:
                strand = -1
            LTR_feature.sub_features.append(SeqFeature(FeatureLocation(transposon['pbs'][0], transposon['pbs'][1]), 
                type='LTR_pbs', strand=1))
        features.append(LTR_feature)        

    #domain type: INT, GAG, AP, RT, RNaseH, Unkwnown
    '''
    domain_features = {}
    domain_features['INT'] = SeqFeature(FeatureLocation(0, len(gene['sequence'])), type='INT', strand=0)
    domain_features['GAG'] = SeqFeature(FeatureLocation(0, len(gene['sequence'])), type='GAG', strand=0)
    domain_features['AP'] = SeqFeature(FeatureLocation(0, len(gene['sequence'])), type='AP', strand=0)
    domain_features['RNaseH'] = SeqFeature(FeatureLocation(0, len(gene['sequence'])), type='RNaseH', strand=0)
    domain_features['RT'] = SeqFeature(FeatureLocation(0, len(gene['sequence'])), type='RT', strand=0)
    domain_features['XX'] = SeqFeature(FeatureLocation(0, len(gene['sequence'])), type='XX', strand=0)
    for df in domain_features:
        domain_features[df].sub_features = []
    #insert domains
    for domain in gene['domains']:
        strand = 1
        x,y = domain['location'][0], domain['location'][1]
        if domain['location'][0] > domain['location'][1]:
            strand = -1
            x,y = y,x
        domain_type = domain['Type'].split('_')[0]
        if domain_type not in ['INT', 'GAG', 'AP', 'RT', 'RNaseH']:
            domain_type = 'XX'
        domain_features[domain_type].sub_features.append(SeqFeature(FeatureLocation(x, y), type=domain_type, strand=strand))

    for df in domain_features:
        features.append(domain_features[df])
    '''
    rec.features = features
    #possible other formats as well
    with open('data/' + gene['id'] + '.gff', 'w+') as out_handle:
        gff = GFF.write([rec], out_handle)

def visualize(gene):
    #geneToGFF(gene)

    #Style
    style = Style()
    style.load_file('gt.style')

    #Feature index
    feature_index = FeatureIndexMemory()

    #add gff file
    feature_index.add_gff3file('data/' + gene['id'] + '.gff')

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
    pngfile = 'data/' + gene['id'] + '.png'
    canvas.to_file(pngfile)
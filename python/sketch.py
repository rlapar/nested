from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from gt.core import *
from gt.extended import *
from gt.annotationsketch import *
from gt.annotationsketch.custom_track import CustomTrack
from gt.core.gtrange import Range
import sys

style = {}

def geneToGFF(gene):
    #insert gene
    rec = SeqRecord(gene['Sequence'], gene['ID'])
    qualifiers = {'ID': gene['ID']}
    features = []
    features.append(SeqFeature(FeatureLocation(0, len(gene['Sequence'])), type='gene', strand=1, qualifiers=qualifiers))
    features[0].sub_features = []

    #insert LTR finder results
    for transposon in gene['LTR_finder']:
        features[0].sub_features.append(SeqFeature(FeatureLocation(transposon['Location'][0], transposon['Location'][1]), 
            type='LTR_finder_TE', strand=1))        

    rec.features = features
    #possible other formats as well
    with open('data/' + gene['ID'] + '.gff', 'w+') as out_handle:
        gff = GFF.write([rec], out_handle)

def visualize(gene):
    geneToGFF(gene)

    #Style
    style = Style()
    style.load_file('gt.style')

    #Feature index
    feature_index = FeatureIndexMemory()

    #add gff file
    feature_index.add_gff3file('data/' + gene['ID'] + '.gff')

    #create diagram for first sequence ID
    seqid = feature_index.get_first_seqid()
    range = feature_index.get_range_for_seqid(seqid)
    diagram = Diagram.from_index(feature_index, seqid, range, style)

    #create layout
    layout = Layout(diagram, 600, style)
    height = layout.get_height()

    #create canvas
    canvas = CanvasCairoFile(style, 600, height)

    #sketch layout on canvas
    layout.sketch(canvas)

    #write to file
    pngfile = 'data/' + gene['ID'] + '.png'
    canvas.to_file(pngfile)
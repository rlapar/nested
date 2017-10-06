from gt.core import *
from gt.extended import *
from gt.annotationsketch import *
from gt.annotationsketch.custom_track import CustomTrack
from gt.core.gtrange import Range
import sys

style = {}

def visualize(transposon):
    print transposon
    gene = FeatureNode.create_new('test', "gene", 0, 60000, "+")
    t = FeatureNode.create_new('test', "exon", transposon['Location'][0], transposon['Location'][1], "+")
    gene.add_child(t)

    pngfile = 'test.png'

    style = Style()

    diagram = Diagram.from_array([gene], Range(1, 1000), style)

    layout = Layout(diagram, 600, style)
    height = layout.get_height()
    canvas = CanvasCairoFile(style, 600, height)
    layout.sketch(canvas)

    canvas.to_file(pngfile)
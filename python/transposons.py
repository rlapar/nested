import math 
import copy
import pprint

from python import config
from python import ltr as pythonLtr
#from python import progressBar
from python import sketch as pythonSketch
from python import domains as pythonBlast
#from python import geneGraph
from python import intervalTree as it

def expandIntervals(intervals):
    for i in range(len(intervals) - 2, -1, -1):
        for j in range(i+1, len(intervals)):
            if intervals[i][0] <= intervals[j][0]:
                intervals[j][0] += (intervals[i][1] - intervals[i][0])
            if intervals[i][0] <= intervals[j][1]:
                intervals[j][1] += (intervals[i][1] - intervals[i][0])
        
    return intervals


def findNestedTranspononsTree(fasta_sequences):
    nested = findNestedTransposons(fasta_sequences)
    genes = pythonLtr.getGeneDict(fasta_sequences)
    for g in genes:
        genes[g]['nested'] = nested[g]
    #reconstruct
    #ITERATE FROM THE BACK AND RECONSTRUCT TREE
    for g in genes:
        genes[g]['nested']['intervals'] = expandIntervals(genes[g]['nested']['intervals'])
        print g, ':', genes[g]['nested']['intervals']
        genes[g]['nested']['tree'] = it.IntervalTree(genes[g]['nested']['intervals'], genes[g]['nested']['transposons']).toDict()
                
    pythonSketch.sketch(genes)  
    

def findNestedTransposons(fasta_sequences):
    #print '______________________'
    genes = pythonLtr.getGeneDict(fasta_sequences)

    #Find LTR pairs using LTR_finder
    genes = pythonLtr.runLtrFinder(genes)

    #Add features (other LTR_finding, BLAST similarity, protein domains)
    genes = pythonBlast.runBlastx(genes)

    #Sketch
    #pythonSketch.sketch(genes)

    #FIND BEST TRANSPOSON
    genes = findBestCandidate(genes)

    #SAVE COORDINATES FOR BACKTRACK
    nested = {}
    for g in genes:
        if math.isnan(genes[g]['best_transposon']['location'][0]):
            nested[g] = {
                'intervals': [],
                'transposons': []
            }
        else:
            nested[g] = {
                'intervals': [genes[g]['best_transposon']['location']],
                'transposons': [genes[g]['best_transposon']]
            }

    #CROP SEQUENCES AND RECURSIVELY CALL
    cropped = False
    cropped_sequences = []
    for sequence in fasta_sequences:
        if not math.isnan(genes[sequence.id]['best_transposon']['location'][0]):
            #print 'Cropped: ' +  str(genes[sequence.id]['best_transposon'])
            cropped = True
            location = genes[sequence.id]['best_transposon']['location']
            cropped_sequences.append(sequence[:(location[0] - 1)] + sequence[(location[1] + 1):])

    #RECURSIVELY CALL
    if cropped:
        nested_cropped = findNestedTransposons(cropped_sequences)
        for a in nested_cropped:
            nested[a]['intervals'] += nested_cropped[a]['intervals']
            nested[a]['transposons'] += nested_cropped[a]['transposons']

    return nested
    

def findBestCandidate(genes):
    for gene in genes:
        for transposon in genes[gene]['ltr_finder']:
            transposon['evaluated'] = evaluateTransposon(transposon, genes[gene]['domains'])

    #get best evaluated transposon
    for gene in genes:
        if not len(genes[gene]['ltr_finder']):
            genes[gene]['best_transposon'] = {
                'location': [float('nan'), float('nan')]
            }
        else:
            best_transposon = max(genes[gene]['ltr_finder'], key=lambda d: d['evaluated']['score'])
            genes[gene]['best_transposon'] = best_transposon
            
    return genes

def evaluateTransposon(transposon, domains):
    score = 0
    features = {
        'domains': [],
        'pbs': None,
        'ppt': None
    }
    #print 'Evaluation transposon: ' + transposon['seqid'] + ' ' + str(transposon['location'])
    
    #Look for domains
    for domain_type in domains:
        for domain in domains[domain_type]:
            domain_location = domain['location']
            if domain_location[0] > domain_location[1]:
                domain_location = [domain_location[1], domain_location[0]]
            #found domain
            if domain_location[0] >= transposon['location'][0] and domain_location[1] <= transposon['location'][1]:
                score += domain['score']
                features['domains'].append(domain)
    
    #Check PBS
    if not math.isnan(transposon['pbs'][0]):
        score += 1000
        features['pbs'] = transposon['pbs']

    #Check PPT
    if not math.isnan(transposon['ppt'][0]):
        score += 1000
        features['ppt'] = transposon['ppt']

    score /= float(transposon['location'][1] - transposon['location'][0])
    
    return {
        'score': score,
        'features': features
    }

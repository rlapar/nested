import os
import subprocess

from Bio import SeqIO
from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from python import config
from python import intervals

def createGFF(id, sequence, nestedList):
    #find closest parent
    parents = [-1] * len(nestedList)
    for i in range(len(parents) - 1):
        for j in range(i + 1, len(parents)):
            if intervals.contains(nestedList[j].location, nestedList[i].location):
                parents[i] = j
                break
    #append direct children
    directChildren = [[] for i in range(len(nestedList))]
    for i in reversed(range(len(directChildren))):
        parent = parents[i]
        if parent != -1:
            directChildren[parent].append(nestedList[i].location)

    #GFF
    rec = SeqRecord(sequence, id)
    features = []

    for i in range(len(nestedList)):
        #insert BaseLine
        features.append(SeqFeature(
            FeatureLocation(nestedList[i].location[0], nestedList[i].location[1]),
            type='te_base',
            strand=0,
            qualifiers={'ID': 'TE_BASE {}'.format(i)}
        ))

        #insert element cropped by its children
        #save element fasta
        subseq = Seq('')
        children = directChildren[i]
        cropped = intervals.crop(nestedList[i].location, children)
        for subinterval in cropped:
            subseq += sequence[subinterval[0]: (subinterval[1] + 1)]
            features.append(SeqFeature(
                FeatureLocation(subinterval[0], subinterval[1]),
                type='te',
                strand=0,
                qualifiers={'ID': 'TE {}'.format(i), 'Parent': 'TE_BASE {}'.format(i)}
            ))

        with open('data/{}/TE/{}.fa'.format(id, i), 'w') as outfile:
            SeqIO.write(
                SeqRecord(subseq, id=id, description='TE-{}'.format(i)),
                outfile,
                'fasta'
            )

        #insert domains
        j = 0
        for domain in nestedList[i].features['domains']:
            sign = (lambda x: x and (1, -1)[x < 0])(domain.frame[0])
            if sign < 0:
                domain.location = [domain.location[1], domain.location[0]]

            overlap = [x for x in children if intervals.contains(domain.location, x)]
            cropped_domain = intervals.crop(domain.location, overlap)
            for part in cropped_domain:
                features.append(SeqFeature(
                    FeatureLocation(part[0], part[1]),
                    type=domain.type,
                    strand=sign,
                    qualifiers={'ID': 'DOMAIN {}-{}'.format(i,j),'Parent': 'TE_BASE {}'.format(i)}
                ))
            j += 1

    rec.features = features
    
    #other formats possible
    with open('data/{}/{}.gff'.format(id, id), 'w+') as out_handle:
        GFF.write([rec], out_handle)

def sketch(id):
    null = open(os.devnull, 'w')
    process = subprocess.check_output([config.gt_sketch_path] + config.gt_sketch_args + ['data/{}/{}.png'.format(id, id)] + ['data/{}/{}.gff'.format(id, id)], stderr=null)
    



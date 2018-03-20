#!/usr/bin/env python3

import os
import subprocess

from Bio import SeqIO
from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from nested.config import config
from nested.utils import intervals

DEFAULT_DIRPATH = 'data'

class Sketcher(object):
    def __init__(self):
        pass

    def _create_dirs(self, element_id, dirpath):
        if not os.path.exists('{}/{}'.format(dirpath, element_id)):
            os.makedirs('{}/{}'.format(dirpath, element_id))
        if not os.path.exists('{}/{}/TE'.format(dirpath, element_id)):
            os.makedirs('{}/{}/TE'.format(dirpath, element_id))

    def create_gff(self, nested_element, dirpath=None):
        if not dirpath:
            dirpath = DEFAULT_DIRPATH

        self._create_dirs(nested_element.id, dirpath)

        #find closes parent
        nl = nested_element.nested_list
        parents = [-1] * len(nl)
        for i in range(len(parents) - 1):
            for j in range(i + 1, len(parents)):
                if intervals.contains(nl[j].location, nl[i].location):
                    parents[i] = j
                    break

        #append direct children
        direct_children = [[] for i in range(len(nl))]
        for i in reversed(range(len(direct_children))):
            parent = parents[i]
            if parent != -1:
                direct_children[parent].append(nl[i].location)

        #GFF
        rec = SeqRecord(nested_element.sequence, nested_element.id)
        features = []

        for i in range(len(nl)):
            #insert baseline
            features.append(SeqFeature(
                FeatureLocation(nl[i].location[0], nl[i].location[1]),
                type='te_base',
                strand=0,
                qualifiers={'ID': 'TE_BASE {}'.format(i)}
            ))

            #insert element cropped by its children
            subseq = Seq('')
            children = direct_children[i]
            cropped = intervals.crop(nl[i].location, children)
            for subinterval in cropped:
                subseq += nested_element.sequence[subinterval[0]: (subinterval[1] + 1)]
                features.append(SeqFeature(
                    FeatureLocation(subinterval[0], subinterval[1]),
                    type='te',
                    strand=0,
                    qualifiers={'ID': 'TE {}'.format(i), 'Parent': 'TE_BASE {}'.format(i)}
                ))

            #save transposon fasta
            with open('{}/{}/TE/{}.fa'.format(dirpath, nested_element.id, i), 'w') as fasta_out:
                SeqIO.write(
                    SeqRecord(subseq, id=nested_element.id, description='TE-{}'.format(i)),
                    fasta_out,
                    'fasta'
                )

            #insert domains
            if 'domains' in nl[i].features:
                j = 0
                for domain in nl[i].features['domains']:
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
                            qualifiers={'ID': 'DOMAIN {}-{}'.format(i, j),'Parent': 'TE_BASE {}'.format(i)}
                        ))
                    j += 1

        #FOR END
        
        rec.features = features

        #create GFF
        gff_filepath = '{}/{}/{}.gff'.format(dirpath, nested_element.id, nested_element.id)
        with open(gff_filepath, 'w+') as gff_out:
            GFF.write([rec], gff_out)

        return gff_filepath

    #def sketch(self, filepath=None):
    def sketch(self, id, dirpath=None):
        if not dirpath:
            dirpath = DEFAULT_DIRPATH
        #self._create_dirs(id, dirpath)
        null = open(os.devnull, 'w')
        process = subprocess.check_output(
            [config.gt_sketch_path] + 
            config.gt_sketch_args + 
            ['{}/{}/{}.png'.format(dirpath, id, id)] + 
            ['{}/{}/{}.gff'.format(dirpath, id, id)], stderr=null)
        

#!/usr/bin/env python3
import os
import math

from Bio import SeqIO
from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from nested.utils import intervals
from nested.config.config import config

format_dict = config['gff_format']

class GFFMaker(object):
    def __init__(self):
        pass

    def _create_dirs(self, element_id, dirpath):
        if not os.path.exists('{}/{}'.format(dirpath, element_id)):
            os.makedirs('{}/{}'.format(dirpath, element_id))
        if not os.path.exists('{}/{}/TE'.format(dirpath, element_id)):
            os.makedirs('{}/{}/TE'.format(dirpath, element_id))

    def create_gff(self,
                   nested_element,
                   dirpath,
                   format='default'):

        if format not in format_dict:
            format = 'default'

        self._create_dirs(nested_element.id, dirpath)

        #TODO move to separate method
        # find closest parent
        nl = nested_element.nested_list
        parents = [-1] * len(nl)
        for i in range(len(parents) - 1):
            for j in range(i + 1, len(parents)):
                if intervals.contains(nl[j].location, nl[i].location):
                    parents[i] = j
                    break

        # append direct children
        direct_children = [[] for i in range(len(nl))]
        for i in reversed(range(len(direct_children))):
            parent = parents[i]
            if parent != -1:
                direct_children[parent].append(nl[i].location)

        # GFF
        rec = SeqRecord(nested_element.sequence, nested_element.id)
        features = []

        for i in range(len(nl)):
            #insert baseline
            base_type = format_dict[format]['te_base'] if format != 'default' else 'te_base'
            features.append(SeqFeature(
                FeatureLocation(
                    nl[i].location[0], nl[i].location[1]),
                    type=base_type,
                    strand=0,
                    qualifiers={
                        'name': 'TE_BASE {}'.format(i),
                        'ID': 'TE_BASE {}'.format(i)
                    }
            ))

            #insert element cropped by its children
            subseq = Seq('')
            children = direct_children[i]
            cropped = intervals.crop(nl[i].location, children)
            for subinterval in cropped:
                subseq += nested_element.sequence[subinterval[0]: (subinterval[1] + 1)]
                te_type = format_dict[format]['te'] if format != 'default' else 'te'
                features.append(SeqFeature(
                    FeatureLocation(subinterval[0], subinterval[1]),
                    type=te_type,
                    strand=0,
                    qualifiers={
                        'ID': 'TE {}'.format(i), 
                        'name': 'TE {}'.format(i), 
                        'Parent': 'TE_BASE {}'.format(i)
                    }
                ))

            # save transposon fasta
            with open('{}/{}/TE/{}.fa'.format(dirpath, nested_element.id, i), 'w') as fasta_out:
                SeqIO.write(
                    SeqRecord(subseq, id='{}|TE-{}'.format(nested_element.id, i), description='Cropped nested retrotransposon'),
                    fasta_out,
                    'fasta'
                )

            # insert domains
            if 'domains' in nl[i].features:
                j = 0
                for domain in nl[i].features['domains']:
                    domain_location = domain.location
                    sign = (lambda x: x and (1, -1)[x < 0])(domain.frame[0])
                    if sign < 0:
                        domain_location = [domain_location[1], domain_location[0]]
                    overlap = [x for x in children if intervals.contains(domain_location, x)]
                    cropped_domain = intervals.crop(domain_location, overlap)
                    for part in cropped_domain:
                        domain_type = format_dict[format]['domain'] if format != 'default' else domain.type
                        features.append(SeqFeature(
                            FeatureLocation(part[0], part[1]),
                            type=domain_type,
                            strand=sign,
                            qualifiers={
                                'ID': 'DOMAIN {}-{}'.format(i, j), 
                                'name': domain.type,
                                'Parent': 'TE_BASE {}'.format(i)
                            }
                        ))
                    j += 1

            #insert pbs,ppt
            if 'pbs' in nl[i].features and not math.isnan(nl[i].features['pbs'][0]):
                pbs_tybe = format_dict[format]['pbs'] if format != 'default' else 'pbs'
                features.append(SeqFeature(
                    FeatureLocation(nl[i].features['pbs'][0], nl[i].features['pbs'][1]),
                    type=pbs_tybe,
                    strand=0,
                    qualifiers={
                        'ID': 'PBS {}'.format(i),
                        'name': 'pbs',
                        'Parent': 'TE_BASE {}'.format(i)
                    }
                ))

            if 'ppt' in nl[i].features and not math.isnan(nl[i].features['ppt'][0]):
                ppt_type = format_dict[format]['ppt'] if format != 'default' else 'ppt'
                features.append(SeqFeature(
                    FeatureLocation(nl[i].features['ppt'][0], nl[i].features['ppt'][1]),
                    type=ppt_type,
                    strand=0,
                    qualifiers={
                        'ID': 'PPT {}'.format(i),
                        'name': 'ppt',
                        'Parent': 'TE_BASE {}'.format(i)
                    }
                ))

            #insert ltrs            
            ltr_type = format_dict[format]['ltr'] if format != 'default' else 'ltr'
            features.append(SeqFeature(
                FeatureLocation(nl[i].ltr_right_location[0], nl[i].ltr_right_location[1]),
                type=ltr_type,
                strand=0,
                qualifiers={
                    'ID': 'LTR RIGHT {}'.format(i),
                    'name': 'ltr right',
                    'Parent': 'TE_BASE {}'.format(i)
                }
            ))
            features.append(SeqFeature(
                FeatureLocation(nl[i].ltr_left_location[0], nl[i].ltr_left_location[1]),
                type=ltr_type,
                strand=0,
                qualifiers={
                    'ID': 'LTR LEFT {}'.format(i),
                    'name': 'ltr left',
                    'Parent': 'TE_BASE {}'.format(i)
                }
            ))

        # FOR END
        rec.features = features

        #create GFF
        filename = '{}/{}/{}'.format(dirpath, nested_element.id, nested_element.id)
        if format != 'default':
            filename += '_{}'.format(format)
        gff_filepath = '{}.gff'.format(filename)
        with open(gff_filepath, 'w+') as gff_out:
            GFF.write([rec], gff_out)

        return gff_filepath

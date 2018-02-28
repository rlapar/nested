#!/usr/bin/env python3

from nested.utils import intervals
from nested.core.gene import Gene
from nested.core.nested_element import NestedElement

class Nester(object):
    """Class represent nesting in sequence, recursivelly find best evaluated transposon and crop it until no new transposons found

    Attributes:
        seqid (str): id of sequence
        sequence (Bio.Seq.Seq): sequence
        nested_element (NestedElement): nested element
    """     
    def __init__(self, sequence):
        self.seqid = sequence.id
        self.sequence = sequence.seq
        self.nested_element = None
        self._find_nesting()

    def _find_nesting(self):
        nested_list = self._get_unexpanded_transposon_list(self.sequence) #find list of nested transposons
        nested_list = self._expand_transposon_list(nested_list) 
        self.nested_element = NestedElement(self.seqid, self.sequence, nested_list)

    def _get_unexpanded_transposon_list(self, sequence): #recursivelly find and crop best evaluated transposon, return unexpanded list of found transposons
        gene = Gene(self.seqid, sequence)
        best_candidate = gene.get_best_candidate()
        if not best_candidate:
            return []
        nested_list = [best_candidate]
        #crop TE and call recursivelly
        cropped_sequence = sequence[:(best_candidate.location[0] - 1)] + sequence[(best_candidate.location[1] + 1):]
        nested_list += self._get_unexpanded_transposon_list(cropped_sequence)
        return nested_list

    def _expand_transposon_list(self, nested_list): #backwards expanding of intervals according to previously found and cropped elements
        for i in reversed(range(len(nested_list) - 1)):
            for j in range(i + 1, len(nested_list)):
                nested_list[j].location = intervals.expand(nested_list[i].location, nested_list[j].location)
                for domain in nested_list[j].features['domains']:
                    domain.location = intervals.expand(nested_list[i].location, domain.location)
                nested_list[j].features['ppt'] = intervals.expand(nested_list[i].location, nested_list[j].features['ppt'])
                nested_list[j].features['pbs'] = intervals.expand(nested_list[i].location, nested_list[j].features['pbs'])
        return nested_list

    

    
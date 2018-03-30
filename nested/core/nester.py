#!/usr/bin/env python3

import os
import re

from nested.utils import intervals
from nested.core.gene import Gene
from nested.core.nested_element import NestedElement
from nested.logging.logger import NesterLogger
from nested.config import config

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
        self._iteration = 0
        self._logger = NesterLogger('{}/{}'.format(config.logdir, self.seqid))
        self._find_nesting()                    

    def _find_nesting(self):
        nested_list = self._get_unexpanded_transposon_list(self.sequence) #find list of nested transposons
        nested_list = self._expand_transposon_list(nested_list) 
        self.nested_element = NestedElement(self.seqid, self.sequence, nested_list)

    def _get_unexpanded_transposon_list(self, sequence): #recursivelly find and crop best evaluated transposon, return unexpanded list of found transposons
        self._iteration += 1
        gene = Gene(self.seqid, sequence)
        candidates = gene.get_candidates_above_threshold(threshold=0)
        if not candidates:
            best_candidate = gene.get_best_candidate()
            if not best_candidate:
                return []
            candidates = [best_candidate]

        #sort by score (from highest)
        candidates.sort(key=lambda x: x.score, reverse=True)        
        #remove intersections
        nested_list = []
        for candidate in candidates:
            choose = True
            for element in nested_list:
                if intervals.intersect(candidate.location, element.location):
                    choose = False
                    break
            for other_candidate in candidates:
                if candidate.location == other_candidate.location:
                    continue
                if intervals.contains(candidate.location, other_candidate.location):
                    choose = False
                    break
            if choose:
                nested_list.append(candidate)
        #sort by location (reverse)
        nested_list.sort(key=lambda x:x.location[0], reverse=True)
        #crop from sequence
        cropped_sequence = sequence

        for element in nested_list:
            cropped_sequence = cropped_sequence[:(element.location[0] - 1)] + cropped_sequence[(element.location[1] + 1):]

        #LOG
        self._logger.log_iteration(self._iteration, nested_list)
        
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

    

    
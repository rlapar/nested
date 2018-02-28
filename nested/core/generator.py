#!/usr/bin/env python3

import random

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from nested.core.te import TE
from nested.core.nested_element import NestedElement
from nested.utils import intervals

class Generator(object):
    """Class used to generate artificial nested elements for testing purposes

    Attributes:
        sourceFile (str): path to source database of TE's in fasta format
    """

    def __init__(self, source_db):
        self._source_db = source_db
        self.elements = []

    def generate_random_nested_element(self, baselength=0, number_of_iterations=1, filter_string='LTR'):
        sequences = list(SeqIO.parse(open(self._source_db), 'fasta'))
        sequences = list(filter(lambda x: filter_string in x.id, sequences))

        #random initial sequence
        element_sequence = ''.join([random.choice('atgc') for x in range(baselength)])
        nested_intervals = [[0, len(element_sequence) - 1]]

        for i in range(number_of_iterations):
            chosen_id = random.randint(0, len(sequences) - 1)
            chosen_seq = sequences[chosen_id].seq
            #print('Chosen sequence ({}) with length {}'.format(sequences[chosen_id].id, len(chosen_seq)))
            insert_position = random.randint(0, len(element_sequence) - 1)
            nested_intervals = [[insert_position, insert_position + len(chosen_seq) - 1]] + nested_intervals
            element_sequence = element_sequence[:insert_position] + chosen_seq + element_sequence[insert_position:]

        #expand intervals
        nested_intervals = intervals.expand_list(nested_intervals)
        nested_tes = []
        for interval in nested_intervals:
            nested_tes.append(TE(
                location=interval
            ))
        
        element_id = 'GENERATED_{}'.format(len(self.elements) + 1)
        element = NestedElement(element_id, element_sequence, nested_tes)
        self.elements.append(element)       

    def generate_random_nested_elements(self, number_of_elements=1, baselength=0, number_of_iterations=1, filter_string='LTR'):
        for i in range(number_of_elements):
            self.generate_random_nested_element(baselength=baselength, number_of_iterations=number_of_iterations, filter_string=filter_string)

    def save_elements_to_fasta(self, filepath):
        sequences = [SeqRecord(x.sequence, id=x.id, description='') for x in self.elements]
        with open(filepath, 'w') as output_handle:
            SeqIO.write(sequences, output_handle, 'fasta')
        
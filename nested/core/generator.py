#!/usr/bin/env python3

import random

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from nested.core.gene import Gene
from nested.core.te import TE
from nested.core.nested_element import NestedElement
from nested.utils import intervals

class Generator(object):
    """Class used to generate artificial nested elements for testing purposes

    Attributes:
        source_db (str): path to source database of TE's in fasta format
    """

    def __init__(self, source_db):
        self._source_db = source_db
        self.elements = []

    def filter_db(self, filtered_db_path, filter_string='LTR', ltr_offset=100, verbose=False):
        sequences = list(SeqIO.parse(open(self._source_db), 'fasta'))
        sequences = list(filter(lambda x: filter_string in x.id, sequences))
        with open(filtered_db_path, 'w+') as filtered_file:
            i = 1
            number_of_entries = len(sequences)
            for sequence in sequences:
                if verbose:
                    print('Filtering entries: {}/{} [{:.2f}%]'.format(i+1, number_of_entries, (100*float(i+1))/number_of_entries), end='\r')
                    i += 1
                sequence.id = sequence.id.replace('/', '--')
                gene = Gene(sequence.id, sequence.seq)
                best_candidate = gene.get_best_candidate()
                if (best_candidate and 
                        best_candidate.location[0] - 1 <= 100 and
                        len(sequence) - best_candidate.location[1] <= 100):                
                    SeqIO.write(sequence, filtered_file, 'fasta')           
        if verbose:
            print('Filtering entries: DONE')

    def generate_random_nested_element(self, baselength=10, number_of_iterations=1):
        sequences = list(SeqIO.parse(open(self._source_db), 'fasta'))

        #random initial sequence
        element_sequence = ''.join([random.choice('atgc') for x in range(baselength)])
        nested_intervals = [[0, len(element_sequence) - 1]]

        for i in range(number_of_iterations):
            chosen_id = random.randint(0, len(sequences) - 1)
            chosen_seq = sequences[chosen_id].seq
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

    def generate_random_nested_elements(self, number_of_elements=1, baselength=10, number_of_iterations=1):
        for i in range(number_of_elements):
            self.generate_random_nested_element(baselength=baselength, number_of_iterations=number_of_iterations)

    def save_elements_to_fasta(self, filepath):
        sequences = [SeqRecord(x.sequence, id=x.id, description='') for x in self.elements]
        with open(filepath, 'w') as output_handle:
            SeqIO.write(sequences, output_handle, 'fasta')
        
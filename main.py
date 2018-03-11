#!/usr/bin/env python3

import sys
import os
import shutil
import click
from datetime import datetime
from subprocess import CalledProcessError

from Bio import SeqIO

from nested.core.nester import Nester
from nested.core.generator import Generator
from nested.output.sketcher import Sketcher

def setup():
	if not os.path.exists('data'):
		os.makedirs('data')

	if not os.path.exists('generated_data'):
		os.makedirs('generated_data')

	if not os.path.exists('tmp'):
		os.makedirs('tmp')

def cleanup():
	shutil.rmtree('tmp/')

@click.command()
@click.option('--input_fasta', '-i', type=str, help='Input fasta file.')
@click.option('--sketch_only', '-s', is_flag=True, help='If true, nesting is not computed. Genes are sketched only from existing gff files.')
@click.option('--generate_from_fasta', '-g', type=str, help='Generate nested elements for source_db.')
@click.option('--baselength', '-gb', type=int, help='Baselength for generated elements.')
@click.option('--number_of_iterations', '-gi', type=int, help='Number of iterations in generating.')
@click.option('--filter_string', '-gf', type=str, help='Filter string in generating.')
@click.option('--number_of_elements', '-ge', type=int, help='Number of generated elements.')
def main(input_fasta, sketch_only, generate_from_fasta, baselength, number_of_iterations, filter_string, number_of_elements):
	number_of_errors = 0
	start_time = datetime.now()
	setup()
	sketcher = Sketcher()

	if input_fasta:
		sequences = list(SeqIO.parse(open(input_fasta), 'fasta'))
		for sequence in sequences:
			sequence.id = sequence.id.replace('/', '--')
			seq_start_time = datetime.now()
			print('Processing {}...'.format(sequence.id), end='\r')
			try:						
				if not sketch_only:
					nester = Nester(sequence)
					sketcher.create_gff(nester.nested_element)
				sketcher.sketch(sequence.id)
				seq_end_time = datetime.now()
				print('Processing {}: done [{}]'.format(sequence.id, seq_end_time - seq_start_time))	
			except KeyboardInterrupt:
				raise
			except CalledProcessError:
				number_of_errors += 1
				print('Processing {}: SUBPROCESS ERROR'.format(sequence.id))		
			except:
				number_of_errors += 1
				print('Processing {}: UNEXPECTED ERROR:'.format(sequence.id), sys.exc_info()[0])		
	elif generate_from_fasta:
		generator = Generator(generate_from_fasta)
		params = {}
		if baselength: params['baselength'] = baselength
		if number_of_iterations: params['number_of_iterations'] = number_of_iterations
		if filter_string: params['filter_string'] = filter_string		
		if number_of_elements: params['number_of_elements'] = number_of_elements
		generator.generate_random_nested_elements(**params)
		generator.save_elements_to_fasta('generated_data/generated.fa')
		for element in generator.elements:
			sketcher.create_gff(element, 'generated_data')
			sketcher.sketch(element.id, 'generated_data')


	cleanup()
	endTime = datetime.now()
	print('Total time: {}'.format(endTime - start_time))
	print('Number of errors: {}'.format(number_of_errors))

if __name__ == '__main__':
	main()

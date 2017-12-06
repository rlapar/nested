#!/usr/bin/env python3

import sys
import os
import shutil
import click

from Bio import SeqIO

from python.nester import Nester
from python import sketch

def setup():
	if not os.path.exists('data'):
		os.makedirs('data')

	if not os.path.exists('tmp'):
		os.makedirs('tmp')

def cleanup():
	shutil.rmtree('tmp/')

@click.command()
@click.option('--input_fasta', '-i', required=True, type=str, help='Input fasta file.')
@click.option('--sketch_only', '-s', is_flag=True)
def main(input_fasta, sketch_only):
	setup()

	#RUN
	sequences = list(SeqIO.parse(open(input_fasta), 'fasta'))
	for sequence in sequences:
		try:
			if not os.path.exists('data/{}'.format(sequence.id)):
				os.makedirs('data/{}'.format(sequence.id))
			if not os.path.exists('data/{}/TE'.format(sequence.id)):
				os.makedirs('data/{}/TE'.format(sequence.id))
			print('Running Nester for {}...'.format(sequence.id))
			
			if not sketch_only:
				nester = Nester(sequence)
				#for a in nester.nestedList:
				#	print(a)
				#	for b in a.features['domains']:
				#		print('--------------------')
				#		print(b)
				#	print('##########################')
				sketch.createGFF(sequence.id, sequence.seq, nester.nestedList)
			sketch.sketch(sequence.id)
		except KeyboardInterrupt:
			raise
		#except:
		#	print('Unexpected error:', sys.exc_info()[0])			

	cleanup()

if __name__ == '__main__':
    main()
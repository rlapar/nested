#!/usr/bin/env python3

import os
import shutil
import click

from Bio import SeqIO

from python.nester import Nester

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
		if not os.path.exists('data/{}'.format(sequence.id)):
			os.makedirs('data/{}'.format(sequence.id))
		if not os.path.exists('data/{}/TE'.format(sequence.id)):
			os.makedirs('data/{}/TE'.format(sequence.id))
		print('Running Nester for {}...'.format(sequence.id))
		nester = Nester(sequence)

		#TODO 
		#createGFF nester.nestedList

		#TODO 
		#sketch only without nester

		#remove
		for a in nester.nestedList:
			print(a.location)
		#remove

	cleanup()

if __name__ == '__main__':
    main()
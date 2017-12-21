#!/usr/bin/env python3

import sys
import os
import shutil
import click
from datetime import datetime
from subprocess import CalledProcessError

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
@click.option('--sketch_only', '-s', is_flag=True, help='If true, nesting is not computed. Genes are sketched only from existing gff files.')
def main(input_fasta, sketch_only):
	numberOfErrors = 0
	startTime = datetime.now()
	setup()

	sequences = list(SeqIO.parse(open(input_fasta), 'fasta'))
	for sequence in sequences:
		seqStartTime = datetime.now()
		print('Processing {}...'.format(sequence.id), end='\r')
		try:
			if not os.path.exists('data/{}'.format(sequence.id)):
				os.makedirs('data/{}'.format(sequence.id))
			if not os.path.exists('data/{}/TE'.format(sequence.id)):
				os.makedirs('data/{}/TE'.format(sequence.id))
			
			if not sketch_only:
				nester = Nester(sequence)
				sketch.createGFF(sequence.id, sequence.seq, nester.nestedList)
			sketch.sketch(sequence.id)
			seqEndTime = datetime.now()
			print('Processing {}: done [{}]'.format(sequence.id, seqEndTime - seqStartTime))	
		except KeyboardInterrupt:
			raise
		except CalledProcessError:
			numberOfErrors += 1
			print('Processing {}: SUBPROCESS ERROR'.format(sequence.id))		
		except:
			numberOfErrors += 1
			print('Processing {}: UNEXPECTED ERROR:'.format(sequence.id), sys.exc_info()[0])		

	cleanup()
	endTime = datetime.now()
	print('Total time: {}'.format(endTime - startTime))
	print('Number of errors: {}'.format(numberOfErrors))

if __name__ == '__main__':
	main()

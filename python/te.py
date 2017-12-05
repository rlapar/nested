#!/usr/bin/env python3

import subprocess
import re

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from python import config

class TE(object):
	def __init__(self, entry=None):
		self.ppt = entry['ppt']
		self.pbs = entry['pbs']
		self.location = entry['location']
		self.features = None
		self.score = None

	def __str__(self):
		lines = ['{{location = {},'.format(self.location),
				 ' ppt = {},'.format(self.ppt),
				 ' pbs = {},'.format(self.pbs),
				 ' score = {}}}'.format(self.score)]
		return '\n'.join(lines)

def runLtrFinder(seqid, sequence):
	transposons = []

	#prepare tmp fasta file for sequence
	with open('tmp/tmp.fa', 'w+') as tmp_file:
		SeqIO.write(SeqRecord(sequence, id=seqid),
					tmp_file,
					'fasta')

	#feed stdin to LTR finder
	process = subprocess.Popen([config.ltr_finder_path] + config.ltr_finder_args + ['tmp/tmp.fa'], 
								stdout=subprocess.PIPE, 
								stderr=subprocess.PIPE)    

	stdout, stderr = process.communicate()
	parsedOutput = parseRawOutput(stdout)
	for entry in parsedOutput:
		transposons.append(TE(entry))

	return transposons

def parseRawOutput(output):
	ltrOutput = None
	entries = output.decode('utf-8').split('>Sequence: ')
	for e in entries:
		entryName = e.split(' ')[0]
		if entryName == 'Program':
			continue
		return parseLtrTable(e.split('\n')[1:])

	return None

def parseLtrTable(rawTable):
	re_interval = re.compile('[0-9N]+[-][0-9N]+')
	#re_int = re.compile('[0-9N]+')
	#re_float = re.compile('(\d+\.\d*)|N')
	tHead = rawTable[0].split('\t')
	transposons = []

	for line in rawTable[1:]:
		transposon = {}
		#TO DO: test properly on real LTR tables
		if len(line.split('\t')) == 1: #newline (end of table) 
			return transposons
		attributes = line.split('\t')
		for i in range(len(attributes)):            
			if re_interval.match(attributes[i]):
				if attributes[i][0] == 'N':
					transposon[str.lower(tHead[i])] = [float('nan'), float('nan')] 
				else:
					transposon[str.lower(tHead[i])] = [int(attributes[i].split('-')[0]), int(attributes[i].split('-')[1])]
			else:
				try:
					transposon[str.lower(tHead[i])] = float(attributes[i])
				except:
					transposon[str.lower(tHead[i])] = attributes[i]

		transposons.append(transposon)
                
	return transposons



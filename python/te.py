#!/usr/bin/env python3

import subprocess
import re

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from python import config

class TE(object):
	"""Class representing TE.

	Attributes:
		ppt (list): location of ppt
		pbs (list): location of pbs
		location (list): position of TE in gene
		ltrLeftLocation (list): location of left ltr
		ltrRightLocation (list): location of right ltr
		features (dict): dictionary of features assigned to TE (i.e. list of domains on the optimal path in graph)
		score (float): evaluated score of TE
	"""
	def __init__(self, entry=None):
		self.ppt = entry['ppt']
		self.pbs = entry['pbs']
		self.location = entry['location']
		self.ltrLeftLocation = [
			entry['location'][0],
			entry['location'][0] + int(entry['ltr len'].split(',')[0])
		]
		self.ltrRightLocation = [
			entry['location'][1] - int(entry['ltr len'].split(',')[1]),
			entry['location'][1]
		]
		self.features = None
		self.score = None

	def __str__(self):
		lines = ['{{location = {},'.format(self.location),
				 ' left ltr = {},'.format(self.ltrLeftLocation),
				 ' right ltr = {},'.format(self.ltrRightLocation),
				 ' ppt = {},'.format(self.ppt),
				 ' pbs = {},'.format(self.pbs),
				 #' features = {},'.format(self.features),
				 ' score = {}}}'.format(self.score)]
		return '\n'.join(lines)

"""Run LTR finder on sequence and return list of transposons

Arguments:
	seqid (str): sequence id
	sequence (Bio.Seq.Seq): sequence

Returns:
	list[TE]: list of found ltr pairs as a TE class
"""
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

"""Parse raw output from LTR finder

Arguments:
	output (str): ltr finder raw output (w -2)

Returns:
	dict: parsed raw output
"""
def parseRawOutput(output):
	ltrOutput = None
	entries = output.decode('utf-8').split('>Sequence: ')
	for e in entries:
		entryName = e.split(' ')[0]
		if entryName == 'Program':
			continue
		return parseLtrTable(e.split('\n')[1:])

	return None

"""Parse raw ltr table

Arguments:
	rawTable (str): raw table

Returns:
	dict: parsed table
"""
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



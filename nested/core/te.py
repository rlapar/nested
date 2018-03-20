#!/usr/bin/env python3

import os
import subprocess
import re

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from nested.config import config

class TE(object):
    """Class representing TE. Every location is in format [from, to].

    Attributes:
        ppt (list): location of ppt
        pbs (list): location of pbs
        location (list): position of TE in gene
        ltr_left_location (list): location of left ltr
        ltr_right_location (list): location of right ltr
        features (dict): dictionary of features assigned to TE (i.e. list of domains on the optimal path in graph)
        score (float): evaluated score of TE
    """
    def __init__(self, ppt=None, pbs=None, location=None, ltr_left_location=None, ltr_right_location=None, features={}, score=None):
        self.ppt = ppt
        self.pbs = pbs
        self.location = location
        self.ltr_left_location = ltr_left_location
        self.ltr_right_location = ltr_right_location
        self.features = features
        self.score = score

    def __str__(self):
        lines = ['{{location = {},'.format(self.location),
				 ' left ltr = {},'.format(self.ltr_left_location),
				 ' right ltr = {},'.format(self.ltr_right_location),
				 ' ppt = {},'.format(self.ppt),
				 ' pbs = {},'.format(self.pbs),
				 #' features = {},'.format(self.features),
				 ' score = {}}}'.format(self.score)]
        return '\n'.join(lines)

def run_ltr_finder(seqid, sequence):
    """Run LTR finder on sequence and return list of transposons

    Arguments:
        seqid (str): sequence id
        sequence (Bio.Seq.Seq): sequence
        tmp_dir (str): Auxiliary existing directory

    Returns:
        list[TE]: list of found ltr pairs as a TE class
    """
    transposons = []

    if not os.path.exists('/tmp/nested'):
    	os.makedirs('/tmp/nested')

    if not os.path.exists('/tmp/nested/ltr'):
        	os.makedirs('/tmp/nested/ltr')

    with open('/tmp/nested/ltr/{}.fa'.format(seqid), 'w+') as tmp_file: #prepare tmp fasta file for ltr finder
        SeqIO.write(SeqRecord(sequence, id=seqid),
                    tmp_file,
                    'fasta')

    #call LTR finder and feed stdin to it
    process = subprocess.Popen([config.ltr_finder_path] + config.ltr_finder_args + ['/tmp/nested/ltr/{}.fa'.format(seqid)], 
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.PIPE) 

    #os.remove('/tmp/nested/ltr/{}.fa'.format(seqid))   

    stdout, stderr = process.communicate()
    parsed_output = parse_raw_output(stdout)
    for entry in parsed_output:
        params = {
			'ppt': entry['ppt'], 
			'pbs': entry['pbs'], 
			'location': entry['location'], 
			'ltr_left_location': entry['ltr_left_location'], 
			'ltr_right_location': entry['ltr_right_location']
		}
        transposons.append(TE(**params))
    
    return transposons

def parse_raw_output(raw_output):
    """Parse raw output from LTR finder

    Arguments:
        raw_output (str): ltr finder raw output (args -w 2)

    Returns:
        dict: parsed raw output
    """
    entries = raw_output.decode('utf-8').split('>Sequence: ')
    for e in entries:
        entry_name = e.split(' ')[0]
        if entry_name == 'Program':
            continue
        return parse_ltr_table(e.split('\n')[1:])
    return {}

def parse_ltr_table(raw_table):
	re_interval = re.compile('[0-9N]+[-][0-9N]+')
	#re_int = re.compile('[0-9N]+')
	#re_float = re.compile('(\d+\.\d*)|N')
	t_head = raw_table[0].split('\t')
	transposons = []

	for line in raw_table[1:]:
		transposon = {}
		#TO DO: test properly on real LTR tables
		if len(line.split('\t')) == 1: #newline (end of table) 
			return transposons
		attributes = line.split('\t')
		for i in range(len(attributes)):            
			if re_interval.match(attributes[i]):
				if attributes[i][0] == 'N':
					transposon[str.lower(t_head[i])] = [float('nan'), float('nan')] 
				else:
					transposon[str.lower(t_head[i])] = [int(attributes[i].split('-')[0]), int(attributes[i].split('-')[1])]
			else:
				try:
					transposon[str.lower(t_head[i])] = float(attributes[i])
				except:
					transposon[str.lower(t_head[i])] = attributes[i]

		transposon['ltr_left_location'] = [
			transposon['location'][0],
			transposon['location'][0] + int(transposon['ltr len'].split(',')[0])
		]
		transposon['ltr_right_location'] = [
			transposon['location'][1] - int(transposon['ltr len'].split(',')[1]),
			transposon['location'][1]
		]
		transposons.append(transposon)

	return transposons



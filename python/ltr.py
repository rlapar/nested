import StringIO
import subprocess
import os
import re
import math
#import pprint

from Bio import SeqIO
from python import config

def getFastaFromFile(input_file):
    return list(SeqIO.parse(open(input_file), 'fasta'))

def getGeneDict(fasta_sequences):
    genes = {}    
    for fasta in fasta_sequences:
        genes[fasta.id] = {
            'sequence': fasta.seq,
            'id': fasta.id
        }    
    return genes

def runLtrFinder(genes):
    #print 'Running LTR finder...'
    #createTmpFastaFile
    with open('tmp/tmp.fa', 'w+') as tmp_file:
        for gene in genes:
            tmp_file.write('>' + genes[gene]['id'] + '\n')
            tmp_file.write(str(genes[gene]['sequence']) + '\n')

    #feedStdInput to LTR_finder
    process = subprocess.Popen([config.ltr_finder_path] + config.ltr_finder_args + ['tmp/tmp.fa'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)    
    stdout, stderr = process.communicate()
    sequences = stdout.split('>Sequence: ')
    transposons = {}

    #parseResults
    for s in sequences:
        seq_name = s.split(' ')[0]
        if seq_name == 'Program':
            continue
        transposons[seq_name] = parseLTRTable(s.split('\n')[1:])
        
    #append to genes
    for ltr in genes:
        if ltr in transposons:
            genes[ltr]['ltr_finder'] = transposons[ltr]  
        else:
            genes[ltr]['ltr_finder'] = [] 

    return genes

def parseLTRTable(rawTable):
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
#open LTR finder
#input take fasta file
import sys
import click
import subprocess
import pprint
import re
import math

from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

from python import config
from python import progressBar
from python import sketch

@click.command()
@click.option('--input_fasta', '-i', required=True, type=str, help='Input fasta file.')
#main needs to be on the top
def main(input_fasta):
    genes = {}
    fasta_sequences = list(SeqIO.parse(open(input_fasta), 'fasta'))
    for fasta in fasta_sequences:
        genes[fasta.id] = {
            'Sequence': fasta.seq
        }
        print fasta.id
    
    print 'Running LTR finder...'
    transposons = runLTR(input_fasta) 
    for t in transposons:
        genes[t]['LTR_finder'] = transposons[t]

    #Add features (other LTR_finding, BLAST similarity, protein domains)
    #Sketch

def runLTR(input_fasta):
    #returning transposons
    process = subprocess.Popen([config.ltr_finder_path] + config.ltr_finder_args + [input_fasta], stdout=subprocess.PIPE, stderr=subprocess.PIPE)    
    stdout, stderr = process.communicate()
    sequences = stdout.split('>Sequence: ')
    transposons = {}
    for s in sequences:
        seq_name = s.split(' ')[0]
        if seq_name == 'Program':
            continue
        transposons[seq_name] = parseLTRTable(s.split('\n')[1:])
        
    return transposons

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
                   transposon[tHead[i]] = [float('nan'), float('nan')] 
                else:
                    transposon[tHead[i]] = [int(attributes[i].split('-')[0]), int(attributes[i].split('-')[1])]
            else:
                try:
                    transposon[tHead[i]] = float(attributes[i])
                except:
                    transposon[tHead[i]] = attributes[i]
                
        transposons.append(transposon)
    return transposons

if __name__ == '__main__':
    main()
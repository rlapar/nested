#open LTR finder
#input take fasta file
import sys
import click
import subprocess
import pprint
import re
import math

from python import config

@click.command()
@click.option('--input_fasta', '-i', required=True, type=str, help='Input fasta file.')

def main(input_fasta):
    process = subprocess.Popen([config.ltr_finder_path] + config.ltr_finder_args + [input_fasta], stdout=subprocess.PIPE, stderr=subprocess.PIPE)    
    stdout, stderr = process.communicate()
    sequences = stdout.split('>Sequence: ')
    transposons = {}
    for s in sequences:
        seq_name = s.split(' ')[0]
        if seq_name == 'Program':
            continue
        transposons[seq_name] = parseLTRTable(s.split('\n')[1:])
        
    pprint.pprint(transposons)    
    

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
import sys
import os
import click
import shutil

from python import ltr as pythonLtr
from python import transposons as pythonTransposons
from python import sketch as pythonSketch

@click.command()
@click.option('--input_fasta', '-i', required=True, type=str, help='Input fasta file.')
#main needs to be on the top
def main(input_fasta):
    fasta_sequences = []    
    try:
        fasta_sequences = pythonLtr.getFastaFromFile(input_fasta)
    except IOError as e:
        print(e)
        sys.exit(1)

    if not os.path.exists('tmp'):
        os.makedirs('tmp')

    nester = pythonTransposons.Nester(fasta_sequences)
    nested_transposons = nester.getNested()
    genes = nester.getGenes()
    '''
    for g in nested_transposons:
        print g, 
        for n in nested_transposons[g]:
            print n['location'],
        print ''
    ''' 
    pythonSketch.sketch(nested_transposons, genes)

    #clear tmp
    shutil.rmtree('tmp/')

if __name__ == '__main__':
    main()
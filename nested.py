import sys
import os
import click
import shutil

from python import ltr as pythonLtr
from python import transposons

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

    transposons.findNestedTranspononsTree(fasta_sequences)

    #clear tmp
    shutil.rmtree('tmp/')

if __name__ == '__main__':
    main()
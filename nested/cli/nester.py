#!/usr/bin/env python3

import sys
import click
from datetime import datetime
from subprocess import CalledProcessError

from Bio import SeqIO

from nested.core.nester import Nester
from nested.output.sketcher import Sketcher


@click.command()
@click.argument('input_fasta', required=True, type=click.Path(exists=True))
@click.option('--sketch', '-s', is_flag=True, default=False, help='Sketch output.')
@click.option('--format', '-f', default='genome_browser', help='Format for GFF.')
@click.option('--output_fasta_offset', '-o', default=0, help='Number of bases around the element included in output fasta files.')
# TODO DATA_FOLDER
def main(input_fasta, sketch, format, output_fasta_offset):
    number_of_errors = 0
    start_time = datetime.now()
    sequences = list(SeqIO.parse(open(input_fasta), 'fasta'))
    sketcher = Sketcher()
    for sequence in sequences:
        sequence.id = sequence.id.replace('/', '--')
        seq_start_time = datetime.now()
        strlen = 15
        print('Processing {a}...'.format(a=sequence.id[:strlen]), end='\r')
        try:
            nester = Nester(sequence)
            sketcher.create_gff(nester.nested_element, output_fasta_offset=output_fasta_offset, format=format)
            if sketch:
                if format != 'default':
                    sketcher.create_gff(nester.nested_element, output_fasta_offset=output_fasta_offset)
                sketcher.sketch(sequence.id)
            seq_end_time = datetime.now()
            print('Processing {a}: DONE [{b}]'.format(a=sequence.id[:strlen], b=seq_end_time - seq_start_time))
        except KeyboardInterrupt:
            raise
        except CalledProcessError:
            number_of_errors += 1
            print('Processing {}: SUBPROCESS ERROR'.format(sequence.id[:strlen]))
        except:
            number_of_errors += 1
            print('Processing {}: UNEXPECTED ERROR:'.format(sequence.id[:strlen]), sys.exc_info()[0])

    end_time = datetime.now()
    print('Total time: {}'.format(end_time - start_time))
    print('Number of errors: {}'.format(number_of_errors))


if __name__ == '__main__':
    main()

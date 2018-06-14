#!/usr/bin/env python3

import click
from datetime import datetime

from nested.core.generator import Generator
from nested.output.sketcher import Sketcher

@click.command()
@click.argument('input_db', required=True, type=click.Path(exists=True))
@click.argument('output_db', required=True, type=str)
@click.option('--baselength', '-l', type=int, help='Baselength for generated elements.')
@click.option('--number_of_iterations', '-i', type=int, help='Number of iterations in generating.')
@click.option('--number_of_elements', '-n', type=int, help='Number of generated elements.')
@click.option('--filter', '-f', is_flag=True, default=False, type=str, help='Filter database and create new one with given output db path.')
@click.option('--filter_string', '-s', type=str, help='Filter entries by given string [ONLY RELEVANT WITH -filter OPTION].')
@click.option('--filter_offset', '-o', type=int, help='LTR offset allowed [ONLY RELEVANT WITH -filter OPTION].')
def main(input_db, output_db, baselength, number_of_iterations, number_of_elements, filter, filter_string, filter_offset):
    #number_of_errors = 0
    start_time = datetime.now()
    generator = Generator(input_db)
    if filter:
        params = {}
        if filter_string: params['filter_string'] = filter_string
        if filter_offset: params['ltr_offset'] = filter_offset
        params['verbose'] = True
        generator.filter_db(output_db, **params)
    else:
        params = {}
        if baselength: params['baselength'] = baselength
        if number_of_iterations: params['number_of_iterations'] = number_of_iterations
        if number_of_elements: params['number_of_elements'] = number_of_elements
        generator.generate_random_nested_elements(**params)
        generator.save_elements_to_fasta('generated_data/{}'.format(output_db))
        sketcher = Sketcher()
        for element in generator.elements:
            sketcher.create_gff(element, 'generated_data')
            sketcher.sketch(element.id, 'generated_data')
    end_time = datetime.now()
    print('Total time: {}'.format(end_time - start_time))
    #print('Number of errors: {}'.format(number_of_errors))


if __name__ == '__main__':
    main()
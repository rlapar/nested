#!/usr/bin/env python3

import os
import subprocess

from nested.config.config import config, args_dict_to_list
from nested.output.gff import GFFMaker

DEFAULT_DIRPATH = 'data'

class Sketcher(object):
    def __init__(self):
        self._gff_maker = GFFMaker()
        self._gff_path = ''

    def create_gff(self, nested_element, dirpath=None, format='default'):
        if not dirpath:
            dirpath = DEFAULT_DIRPATH
        self._gff_maker.create_gff(nested_element, dirpath, format)

    #def sketch(self, filepath=None):
    def sketch(self, id, dirpath=None):
        if not dirpath:
            dirpath = DEFAULT_DIRPATH
        #self._create_dirs(id, dirpath)
        null = open(os.devnull, 'w')
        process = subprocess.check_output(
            [config['gt']['path'], config['gt']['command']] + 
            args_dict_to_list(config['gt']['args']) + 
            ['{}/{}/{}.png'.format(dirpath, id, id)] + 
            ['{}/{}/{}.gff'.format(dirpath, id, id)], stderr=null)

#!/usr/bin/env python3

import os
from datetime import datetime

class Logger(object):
    def __init__(self, destination):
        self._destination = destination
        self._make_dirs()

    def _make_dirs(self):
        dirs = self._destination.split('/')
        for i in range(len(dirs)):
            if not dirs[i]:
                continue
            subpath = '/'.join(dirs[:(i+1)])
            if not os.path.exists(subpath):
        	    os.makedirs(subpath)

    def create_log(self, filepath):
        open('{}/{}'.format(self._destination, filepath), 'w+').close()

    def log_entry(self, filepath, entry):
        with open('{}/{}'.format(self._destination, filepath), 'a') as logfile:
            logfile.write('[{}]\n'.format(datetime.now()))
            logfile.write('{}\n'.format(entry))

class NesterLogger(Logger):
    def __init__(self, destination):
        super().__init__(destination)

    def log_iteration(self, iteration, nested_list):        
        main_logfile_path = 'iterations.log'
        iteration_logfile_path = 'iter_{:02d}.log'.format(iteration)
        if iteration == 1:
            super().create_log(main_logfile_path)
        main_entry = 'ITERATION {}\n'.format(iteration)
        main_entry += '- Cropping {} elements\n'.format(len(nested_list))
        super().log_entry(main_logfile_path, main_entry)

        super().create_log(iteration_logfile_path)
        iteration_entry = ''
        for i in range(len(nested_list)):
            iteration_entry += ('ELEMENT {}\n\n'.format(i+1))
            iteration_entry += ('{}\n\n'.format(nested_list[i]))
            iteration_entry += ('- domains:\n')
            for domain in nested_list[i].features['domains']:
                iteration_entry += ('{}\n'.format(domain)) 
            iteration_entry += '\n\n'
            super().log_entry(iteration_logfile_path, iteration_entry)

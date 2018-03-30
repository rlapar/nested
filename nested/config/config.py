ltr_finder_path = '/home/radovan/lexa/LTR_finder/ltr_finder'
ltr_finder_prosite_path = '/home/radovan/lexa/prosite'
ltr_finder_tRNAdb_path = '/home/radovan/lexa/LTR_finder/tRNAdb/Athal-tRNAs.fa'

ltr_finder_args = ['-l', '32']
ltr_finder_args += ['-S', '3']
ltr_finder_args += ['-o', '2']
ltr_finder_args += ['-p', '15']
ltr_finder_args += ['-g', '80']
ltr_finder_args += ['-G', '5']
ltr_finder_args += ['-T', '2']
ltr_finder_args += ['-w', '2'] #don't change, parsing format required!
ltr_finder_args += ['-a', ltr_finder_prosite_path]
ltr_finder_args += ['-s', ltr_finder_tRNAdb_path]

gt_sketch_path = 'gt'
gt_sketch_args = ['sketch']
gt_sketch_args += ['-width', '1400']
gt_sketch_args += ['-style', '/home/radovan/git/nested/nested/config/gt.style']

blastx_args = {
    'db': '/home/radovan/lexa/blastDB/data/gydb_proteins.fa',
    'outfmt': 5, #don't change, parsing format required!
    'num_threads': 3,
    'dbsize': 90000,
    'word_size': 2,
    'evalue': 1
}

logdir = '/tmp/nested/logs'
#ADD BLAST PARAMETERS

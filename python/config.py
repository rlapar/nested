ltr_finder_path = "/home/radovan/lexa/LTR_finder/ltr_finder"
ltr_finder_prosite_path = "/home/radovan/lexa/prosite"
ltr_finder_tRNAdb_path = "/home/radovan/lexa/LTR_finder/tRNAdb/Athal-tRNAs.fa"

ltr_finder_args = ['-l', '32']
ltr_finder_args += ['-S', '3']
ltr_finder_args += ['-o', '2']
ltr_finder_args += ['-p', '15']
ltr_finder_args += ['-g', '80']
ltr_finder_args += ['-G', '5']
ltr_finder_args += ['-T', '2']
ltr_finder_args += ['-w', '2']
ltr_finder_args += ['-a', ltr_finder_prosite_path]
ltr_finder_args += ['-s', ltr_finder_tRNAdb_path]
 
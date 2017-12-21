# nested

#### External tools
Program is using other external tools. Before usage make sure you have the following tools installed:

- LTR finder (v 1.05+) -
http://www.mybiosoftware.com/ltr_finder-1-0-5-find-full-length-ltr-retrotranspsons-genome-sequences.html
- blastx (v 2.2.31+) -
https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download
- genometools (v 1.5.8) with gt sketch -
https://github.com/genometools/genometools

#### Python3 libraries
Program is running only on Python 3 version and uses non-standard libraries:

- BioPython -
http://biopython.org/
- BCBio GFF -
https://github.com/chapmanb/bcbb/tree/master/gff
- networkx -
https://networkx.github.io/

#### Setting up the config
Before usage, **python/config.py** needs to be set up:

- ltr paths - executable, prosite and tRNAdb
- ltr arguments 
- gt sketch path
- gt sketch arguments
- blastx protein database (our database available at ... )
- blastx arguments

```
Usage: python3 main.py [OPTIONS]

Options:
  -i, --input_fasta TEXT  Input fasta file.  [required]
  -s, --sketch_only       If true, nesting is not computed. Genes are sketched
                          only from existing gff files.
  --help                  Show this message and exit.
```

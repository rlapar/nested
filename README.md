# nested

#### External tools
Program is using other external tools. Before usage make sure you have the following tools installed:

- LTR finder (v 1.05+) -
http://www.mybiosoftware.com/ltr_finder-1-0-5-find-full-length-ltr-retrotranspsons-genome-sequences.html
- blastx (v 2.2.31+) -
https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download
- genometools (v 1.5.8) with gt sketch -
https://github.com/genometools/genometools

#### Install

```
chmod +x setup.sh
sudo ./setup.sh
```

#### Setting up the config
Before usage, **config.yml** needs to be set up:

- ltr paths - executable, prosite and tRNAdb
- ltr arguments 
- gt sketch path
- gt sketch arguments
- blastx protein database (our database available at ... )
- blastx arguments

**Note**: `config.yml` is located by default in directory `/etc/nested/.`

#### Usage
```
Usage: nested-nester [OPTIONS] INPUT_FASTA

Options:
  -s, --sketch_only       If true, nesting is not computed. Genes are sketched
                          only from existing gff files.
  -d, --data_folder TEXT  Output data folder.
  --help                  Show this message and exit.
```

```
Usage: nested-generator [OPTIONS] INPUT_DB OUTPUT_DB

Options:
  -l, --baselength INTEGER        Baselength for generated elements.
  -i, --number_of_iterations INTEGER
                                  Number of iterations in generating.
  -n, --number_of_elements INTEGER
                                  Number of generated elements.
  -f, --filter                    Filter database and create new one with
                                  given output db path.
  -s, --filter_string TEXT        Filter entries by given string [ONLY
                                  RELEVANT WITH -filter OPTION].
  -o, --filter_offset INTEGER     LTR offset allowed [ONLY RELEVANT WITH
                                  -filter OPTION].
  -d, --data_folder TEXT          Output data folder.
  --help                          Show this message and exit.
```

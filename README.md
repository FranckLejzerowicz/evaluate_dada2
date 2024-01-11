# evaluate_dada2

## Description

This denoizing using DADA2 is an exploration of different trim lengths. It 
includes a mock community, clustered at different perc_identity to allow for 
estimates at different levels of uncertainty with respect to the mock. The 
results are written and used for the making of diagnostic plots, contained 
in one single, annotated PDF file.

## Installation
```
pip install git+https://github.com/FranckLejzerowicz/evaluate_dada2.git
```

### Usage

```
Usage: run_dada2 [OPTIONS]

Options:
  -i, --i-fastq-dir TEXT          Folder containing the fastq files
  -m, --i-metadata TEXT           Metadata file (tab-separated) column
                                  `sample_name` match fastq files
  -mi, --i-mock-dir TEXT          Folder containing the mocks sequnces and
                                  clusters
  -mt, --i-mock-tax TEXT          Name of the taxonomy .tsv file in the
                                  `--i-mock-dir` folder
  -r, --p-ranks TEXT              Taxonomy ranks for the mock evaluation
  -c, --p-meta-vars TEXT          Metadata variables to use for checking
                                  samples vs mock
  -r, --p-trim-range INTEGER...   Min, Max, Step for the trimming length
  -l, --p-trim-lengths INTEGER    Trimming lengths
  -lf, --p-f-trim-lengths INTEGER
                                  Trimming lengths (forward reads)
  -lr, --p-r-trim-lengths INTEGER
                                  Trimming lengths (reverse reads)
  -n, --p-n-cores INTEGER         Number of cores for multiprocessing
  --version                       Show the version and exit.
  --help                          Show this message and exit.
```


### Bug Reports

contact `franck.lejzerowicz@gmail.com`
# pycrispr
 *a tool for designing crispr libraries*

 [![Build Status](
 https://travis-ci.org/alaindomissy/pycrispr.svg?branch=master)](https://travis-ci.org/alaindomissy/pycrispr)
 [![Coverage Status](https://coveralls.io/repos/github/alaindomissy/pycrispr/badge.svg?branch=master)](https://coveralls.io/github/alaindomissy/pycrispr?branch=master)
 [![Python 2.7 Status](https://img.shields.io/badge/Python-2.7-brightgreen.svg)](https://img.shields.io/badge/Python-2.7-blue.svg)
 [![Python 3.3 Status](https://img.shields.io/badge/Python-3.3-brightgreen.svg)](https://img.shields.io/badge/Python-3.3-blue.svg)
 [![Python 3.4 Status](https://img.shields.io/badge/Python-3.4-brightgreen.svg)](https://img.shields.io/badge/Python-3.4-blue.svg)
 [![Python 3.5 Status](https://img.shields.io/badge/Python-3.5-brightgreen.svg)](https://img.shields.io/badge/Python-3.5-blue.svg)
  
  
## Installation

```
// sudo apt-get install ncbi-blast+
// results in an outdated version on ubuntu 14.04

        // $ blastn -version
        // blastn: 2.2.28+
        // Package: blast 2.2.28, build Jun  3 2013 11:17:14
// do this instead -actually no, not working either
// conda install -c https://conda.anaconda.org/bioconda blast

// so do this instead
$ wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.3.0+-x64-linux.tar.gz
$ tar xvfp ncbi-blast-2.3.0+-x64-linux.tar.gz
$ export PATH=”$PATH:$HOME/ncbi-blast-2.2.29+/bin”
$ mkdir ./ncbi-blast-2.2.29+/db
$export BLASTDB=”$HOME/ncbi-blast-2.2.29+/db”

        // $ blastn -version
        // blastn: 2.2.31+
        // Package: blast 2.2.31, build Dec  3 2015 17:28:17

sudo apt-get install tree
sudo apt-get install jq

sudo apt-get install bedtools
sudo apt-get install tabix
sudo apt-get install igv 
pip install git+https://github.com/alaindomissy/buffet.git#egg=pycrispr
pip install primer3-py
```

pycrispr runs on Python 2.7, 3.3 and 3.4
 
## Usage

Given a genomic interval, and a rference genome, obtain all candidate crispr guides, 
resulting from enzymatic digestion via cripsr-eating protocol
 
```
enzymatic_protospacers(
    '~/scaffolds_directory/',
    'chr6:136640001-136680000',
    'mm8.fasta'
    )
```

Score candidates, define and rank by yield of good guides all possible clusters of consecutive good guides
 
```
high_specificity_clusters(
    scaffolds_directory,
    genomic_coord, reference,
    genome
    )
```

Advanced options:
 
```
minimum_specificity_clusters(
    scaffolds_directory, genomic_coord, reference, genome,
    chunk_size=25, 
    max_hsps=25,
    ref_substrate_id='chr6',
    low=75, 
    high=75, 
    load_genome=False, 
    howmany=None
    )
```
 
Given a required number of guides, design amplification primers for enough top-yielding good-guides clusters

```
primer_design(scaffolds_directory, genomic_coord, reference, genome, required_number_of_guides)
```
 
## Documentation


## Online access


## Docker image
 
 
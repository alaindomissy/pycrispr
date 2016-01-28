# pycrispr
 *a tool for designing crispr libraries*

 [![Build Status](https://travis-ci.org/alaindomissy/pycrispr.svg?branch=master)](https://travis-ci.org/alaindomissy/pycrispr)
 [![Coverage Status](https://coveralls.io/repos/github/alaindomissy/pycrispr/badge.svg?branch=master)](https://coveralls.io/github/alaindomissy/pycrispr?branch=master)
 [![Python 2.7 Status](https://img.shields.io/badge/Python-2.7-brightgreen.svg)](https://img.shields.io/badge/Python-2.7-blue.svg)
 
## Installation
 

 ```
 pip install git+https://github.com/alaindomissy/buffet.git#egg=pycrispr
 
 pip install primer3-py
 ```

 pycrispr runs on Python 2.7
 
## Usage

Given a genomic interval, and a rference genome, obtain all candidate crispr guides, 
resulting from enzymatic digestion via cripsr-eating protocol
 
```
enzymatic_protospacers('~/scaffolds_directory/', 'chr6:136640001-136680000', 'mm8.fasta')
```

Score candidates, define and rank by yield of good guides all possible clusters of consecutive good guides
 
```
high_specificity_clusters(scaffolds_directory, genomic_coord, reference, genome)
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
 
 
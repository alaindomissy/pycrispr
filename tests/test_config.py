import sys,os
sys.path.insert(0,os.path.abspath(__file__+"/../.."))

from crispr.config import genomes_path, protosp_path, sorted_unique_firstdotsplit, genomes_list, protosp_list


def test_genomes_path():
    assert(genomes_path('mm8') == '/GENOMES/mm8.fasta')


def test_protosp_path():
    assert(protosp_path('mm8') == '/PROTOSP/mm8.prsp.bed')


def test_sorted_unique_firstdotsplit():
    filenames = ['phix.fasta', 'phix.fasta.fai', 'hg38.fasta.fai', 'hg38.fasta','saccer3.fasta', 'saccer3.fasta.fai']
    assert(sorted_unique_firstdotsplit(filenames) == ['hg38', 'phix', 'saccer3'])


def test_genomes_path():
    assert(genomes_list() == ['ecoli', 'hg18', 'hg38', 'lambda', 'mycotube', 'phix', 'saccer3'])


def test_protosp_list():
    assert(protosp_list() == ['ecoli', 'hg18', 'lambda', 'mm8', 'mycotube', 'phix', 'saccer3'])

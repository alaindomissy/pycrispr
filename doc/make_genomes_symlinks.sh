#!/bin/bash

# GENOMES = {
#     'hg38'   : '/genomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa',
#     'hg19'   : '/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa',
#     'hg18'   : '/genomes/Homo_sapiens/UCSC/hg18/Sequence/WholeGenomeFasta/genome.fa',
#     'mm10'   : '/genomes/Homo_sapiens/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa',
#     'mm9'    : '/genomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa',
#     'tair10' : '/genomes/Arabidopsis_thaliana/NCBI/TAIR10/Sequence/WholeGenomeFasta/genome.fa',
#     'saccer3': '/genomes/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/WholeGenomeFasta/genome.fa',
#     'ecoli'  : '/genomes/Escherichia_coli_K_12_DH10B/NCBI/2008-03-17/Sequence/WholeGenomeFasta/genome.fa',
#     'phix'   : '/genomes/PhiX/Illumina/RTA/Sequence/WholeGenomeFasta/genome.fa'
# }


ln -s ~/genomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa ~/GENOMES/hg38.fasta
ln -s ~/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa ~/GENOMES/hg19.fasta
ln -s ~/genomes/Homo_sapiens/UCSC/hg18/Sequence/WholeGenomeFasta/genome.fa ~/GENOMES/hg18.fasta
ln -s ~/genomes/Homo_sapiens/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa ~/GENOMES/mm10.fasta
ln -s ~/genomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa  ~/GENOMES/mm9.fasta
ln -s ~/genomes/Arabidopsis_thaliana/NCBI/TAIR10/Sequence/WholeGenomeFasta/genome.fa ~/GENOMES/tair10.fasta
ln -s ~/genomes/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/WholeGenomeFasta/genome.fa ~/GENOMES/saccer3.fasta
ln -s ~/genomes/Escherichia_coli_K_12_DH10B/NCBI/2008-03-17/Sequence/WholeGenomeFasta/genome.fa ~/GENOMES/ecoli.fasta
ln -s ~/genomes/PhiX/Illumina/RTA/Sequence/WholeGenomeFasta/genome.fa ~/GENOMES/phix.fasta
ls

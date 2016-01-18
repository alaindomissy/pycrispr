import os

from buffet.cut import cut_file


BLASTDB =  os.environ.get('BLASTDB','/genomes/blastdb/')

SCAFFOLDS =  os.environ.get('SCAFFOLDS','/data/scaffolds/mm8/')



datadirpath = "/data/"
#
# main.blast(datadirpath, 'chr6:136640001-136641000_1000.prsp', '/BLASTDB/', 'mm8')
#
#
# from buffet import cut
# DIRECT = '/data/scaffolds/mm8/'
# cut.create_cutbedlines_from_seq_file(DIRECT, 'chr6-1-400', 'fasta')
#
#
# from buffet import digest
# DIRECT = '/data/scaffolds/hg38/'
# digest.digest_file(DIRECT, 'hg38', 'fasta')
#
#
# digest.digest_coord('./', 'chr21:42774475-42905100','mm8', 'fasta')
#


def payload(params_value):

    fastafile = params_value.get('Input.genomic-coord', 'chr6-1-400.fasta')
    # genome = args_value['Input.genome-id']
    # chunk_size = args_value['Input.blast_chunk_size']
    # max_hsps = args_value['Input.blast_max_hsps']

    # args_value['Input.binding_interference_spacing']

    cuts = cut_file(SCAFFOLDS + fastafile)

    return cuts

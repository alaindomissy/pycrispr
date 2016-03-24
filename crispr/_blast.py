from __future__ import absolute_import, division, print_function             # , unicode_literals
import argparse
from .blast import blast
from .config import BLAST_LOG_ON

parser = argparse.ArgumentParser(description='blast a protospacers bed file to shortlist psosible off-targets locations')

parser.add_argument('filename',
                     help='prsp_file_wo_fasta_ext: protospacers file without .fasta extension to use as blast queries')

parser.add_argument('genome',
                     help="genome: name of the blast database to blast the protospacers against")

parser.add_argument('max_hsps',
                     type=int,
                     help="max_hsps: per scaffold hit, maximum number of hsps to return")

parser.add_argument('chunk_size',
                     type=int,
                     help="chunk_size: number of protospacers to blast in each chunk"
                   )

###

vq_arg_group = parser.add_mutually_exclusive_group()

vq_arg_group.add_argument("-v", "--verbose", help="increase output verbosity",
                    action="store_true")

vq_arg_group.add_argument("-q", "--quiet", help="no logs",
                    action="store_true")

###

args = parser.parse_args()
if args.verbose:
    print("verbosity turned on")
if args.quiet:
    print("logs turned off")
    BLAST_LOG_ON = False

nbr_of_chunks = blast(args.filename, args.genome, './', args.max_hsps, args.chunk_size)
print("number of chunks : ", nbr_of_chunks)

###

"""
$ python -m crispr._blast chr6_136640000-136644000_4000 mm8 10 10
"""

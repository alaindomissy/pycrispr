from __future__ import absolute_import, division, print_function    # , unicode_literals
import argparse
from .score import score
from .config import SCORE_LOG_ON

parser = argparse.ArgumentParser()

parser.add_argument('nbrofchunks',
                    type=int,
                     help="nbrofchunks: number of chunks")

parser.add_argument('filename',
                     help='prsp_file_wo_fasta_ext: protospacers fasta file without the .fasta extension' )

parser.add_argument('genome',
                     help="genome: name of the blast database to blast the protospacers against")

parser.add_argument('chunk_size',
                    type=int,
                     help="chunk_size: size of chunks")

parser.add_argument('low',
                    type=int,
                    help="low: ")

parser.add_argument('high',
                    type=int,
                    help="high:")

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
    SCORE_LOG_ON = False

guides = score(args.nbrofchunks, args.filename, args.genome, './', args.chunk_size, False, args.low, args.high)

print('\nPRINT SORTED SCORED GUIDES')
# print(guides.format('fasta'))  # TODO why not this ?
for guide in guides:
    print(guide.format('fasta'))
    # print(guide.format('gb'))   # TODO why this? Locus identifier 'chr:101153-101173(+)' is too long

###

"""
$ python -m crispr._score 5 chr6_136640000-136644000_4000 mm8 10 75 94
"""

from __future__ import absolute_import, division, print_function    # , unicode_literals
import argparse
from .main import digest_blast_score_cluster_prime
# from .config import PRIME_LOG_ON

parser = argparse.ArgumentParser(description='prime etc...')


parser.add_argument('coord',
                     help='cord: scaffold:start-end or scaffold:start_length')

# parser.add_argument('filename',
#                      help='prsp_file_wo_fasta_ext: protospacers fasta file without the .fasta extension' )

parser.add_argument('genome',
                     help="genome: reference genome")

parser.add_argument('max_hsps',
                     type=int,
                     help="max_hsps: per scaffold hit, maximum number of hsps to return")

parser.add_argument('chunk_size',
                     type=int,
                     help="chunk_size: number of protospacers to blast in each chunk"
                   )
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

# if args.verbose:
#     print("verbosity turned on")
# if args.quiet:
#     print("logs turned off")
#     PRIME_LOG_ON = False

digest_blast_score_cluster_prime(args.coord, args.genome, './', args.max_hsps, args.chunk_size, False, args.low, args.high)

###

"""
$ python -m crispr._main chr:101001-109000 ecoli 10 10 75 94

$ python -m crispr._main chr:102000-106999 ecoli 10 10 75 94

"""

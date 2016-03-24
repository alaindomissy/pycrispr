from __future__ import absolute_import, division, print_function    # , unicode_literals
import argparse
from .amplicon import amplicon
# from .config import STRETCH_LOG_ON

parser = argparse.ArgumentParser(description='amplicon etc...')

parser.add_argument('filename',
                     help='prsp_file_wo_fasta_ext: protospacers fasta file without the .fasta extension' )

parser.add_argument('genome',
                     help="genome: name of the blast database to blast the protospacers against")

parser.add_argument('threshold',
                    type=int,
                    help="threshold: ")

###

# vq_arg_group = parser.add_mutually_exclusive_group()
#
# vq_arg_group.add_argument("-v", "--verbose", help="increase output verbosity",
#                     action="store_true")
#
# vq_arg_group.add_argument("-q", "--quiet", help="no logs",
#                     action="store_true")

###

args = parser.parse_args()

# if args.verbose:
#     print("verbosity turned on")
# if args.quiet:
#     print("logs turned off")
#     STRETCH_LOG_ON = False


sorted_amplicons = amplicon(args.filename, './', args.genome, args.threshold)

###

"""
$ python -m crispr._amplicon chr6_136640000-136644000_4000 mm8 94
"""

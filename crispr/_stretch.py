from __future__ import absolute_import, division, print_function    # , unicode_literals
import argparse
from .stretch import stretch, run
from .config import STRETCH_LOG_ON

parser = argparse.ArgumentParser(description='cluster etc...')


parser.add_argument('filename',
                     help='prsp_file_wo_fasta_ext: protospacers fasta file without the .fasta extension' )

parser.add_argument('low',
                    type=int,
                    help="low: ")

parser.add_argument('high',
                    type=int,
                    help="high:")

parser.add_argument('howmany',
                    type=int,
                    help="howmany:")

parser.add_argument('threshold',
                    type=int,
                    help="threshold: ")

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
    STRETCH_LOG_ON = False

guides, stretches = stretch(args.filename, './', args.low, args.high, args.howmany)
#
# print('\nPRINT STRETCHES')
# print(("substrate, start, end, length, good guides, good guides per thousand nucl."))
# for idx, stretch in enumerate(stretches):
#     print(idx, stretch)


guides, runs = run(args.filename, './', args.threshold)

###
"""
$ python -m crispr._stretch chr6_136640000-136644000_4000 75 94 3 94
"""

from __future__ import absolute_import, division, print_function             # , unicode_literals
import argparse
from .digest import digest_genome
from .config import DIGEST_LOG_ON

parser = argparse.ArgumentParser(description='digest a genome')

parser.add_argument('genome',
                     help="genome: reference genome")


vq_arg_group = parser.add_mutually_exclusive_group()

vq_arg_group.add_argument("-v", "--verbose", help="increase output verbosity",
                    action="store_true")

vq_arg_group.add_argument("-q", "--quiet", help="no logs",
                    action="store_true")

args = parser.parse_args()
if args.verbose:
    print("verbosity turned on")
if args.quiet:
    print("logs turned off")
    DIGEST_LOG_ON = False
all()
digest_genome(args.genome)

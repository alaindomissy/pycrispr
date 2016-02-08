from __future__ import absolute_import, division, print_function
# from __future__ import unicode_literals
import argparse
from .digest import digest_genome
from .config import DIGEST_LOG_ON

parser = argparse.ArgumentParser(description='blast a protospacers bed file to shortlist psoiible off-targets locations')

parser.add_argument('genome',
                     help='genome: hg18, hg19. hg38, mm8, mm9, mm10, etc...')

group = parser.add_mutually_exclusive_group()
group.add_argument("-v", "--verbose",
                   help="increase output verbosity",
                    action="store_true")
group.add_argument("-q", "--quiet",
                   help="no logs",
                    action="store_true")

args = parser.parse_args()
if args.verbose:
    print("verbosity turned on")
if args.quiet:
    print("logs turned off")
    DIGEST_LOG_ON = False

digest_genome(args.genome)
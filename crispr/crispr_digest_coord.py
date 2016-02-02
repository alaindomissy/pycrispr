from __future__ import print_function
import argparse
from digest import digest_coord
from config import DIGEST_LOG_ON

parser = argparse.ArgumentParser(description='blast a protospacers bed file to shortlist psoiible off-targets locations')

parser.add_argument('coord',
                     help='cord: scaffold:start-end or scaffold:start_length')

parser.add_argument('genome',
                     help="genome: reference genome")

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

digest_coord(args.coord, args.reference, './')
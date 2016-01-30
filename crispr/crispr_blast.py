from __future__ import print_function
try:
    from io import StringIO        # python3
except ImportError:
    from StringIO import StringIO  # python2
import argparse
from digest import digest_coord

parser = argparse.ArgumentParser()

parser.add_argument('coord',
                     help='scaffold:start-end or scaffold:start_length')

parser.add_argument('reference',
                     help=" reference fasta file")

parser.add_argument("-v", "--verbose", help="increase output verbosity",
                    action="store_true")

parser.add_argument("-q", "--quiet", help="no logs",
                    action="store_true")

args = parser.parse_args()
if args.verbose:
    print("verbosity turned on")
if args.quiet:
    print("logs turned off")
    DIGEST_LOG_ON = False

digest_coord('./', args.coord, args.reference)

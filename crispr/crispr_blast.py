from __future__ import print_function
import argparse
from blast import blast
from config import BLAST_LOG_ON

parser = argparse.ArgumentParser()

parser.add_argument('prsp_file_wo_fasta_ext',
                     help='prsp_file_wo_fasta_ext: protospacers file without .fasta extension to use as blast queries')

parser.add_argument('blastdb',
                     help="blastdb: name of the blast database to blast the protospacers againt")

parser.add_argument('chunk_size',
                     type=int,
                     help="chunk_size: number of protospacers to blast in each chunk"
                   )

parser.add_argument('max_hsps',
                     type=int,
                     help="max_hsps: per scaffold hit, maximum number of hsps to return")

parser.add_argument("-v", "--verbose", help="increase output verbosity",
                    action="store_true")

parser.add_argument("-q", "--quiet", help="no logs",
                    action="store_true")

args = parser.parse_args()
if args.verbose:
    print("verbosity turned on")
if args.quiet:
    print("logs turned off")
    BLAST_LOG_ON = False

blast(args.prsp_file_wo_fasta_ext, args.blastdb, './', args.max_hsps, args.chunk_size)

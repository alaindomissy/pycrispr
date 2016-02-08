from __future__ import absolute_import, division, print_function
# from __future__ import unicode_literals
import argparse
from .cluster import cluster
from .config import CLUSTER_LOG_ON

parser = argparse.ArgumentParser(description='blast a protospacers bed file to shortlist psoiible off-targets locations')


parser.add_argument('blastdb',
                     help="blastdb: name of the blast database to blast the protospacers againt")

parser.add_argument('low',
                    type=int,
                    help="low: ")

parser.add_argument('high',
                    type=int,
                    help="high:")

parser.add_argument('howmany',
                    type=int,
                    help="howmany:")

parser.add_argument('prsp_file_wo_fasta_ext',
                     help='prsp_file_wo_fasta_ext: protospacers fasta file without the .fasta extension' )




parser.add_argument("-v", "--verbose", help="increase output verbosity",
                    action="store_true")

parser.add_argument("-q", "--quiet", help="no logs",
                    action="store_true")

args = parser.parse_args()
if args.verbose:
    print("verbosity turned on")
if args.quiet:
    print("logs turned off")
    CLUSTER_LOG_ON = False



#cluster(guides, './', None, args.low, args.high, args.howmany,args.prsp_file_wo_fasta_ext)

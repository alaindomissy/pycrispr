# import sys,os
# sys.path.insert(0,os.path.abspath(__file__+"/../.."))

from __future__ import absolute_import, division, print_function
# from __future__ import unicode_literals
import argparse
from .config import genomes_list, blastdb_list, protosp_list

parser = argparse.ArgumentParser(description='lists all available reference: genomes or blastdb or protosp ')

parser.add_argument('filetype',
                     help='filetype: genomes or blastdb or protosp' )

# not really used yet
parser.add_argument('action',
                     help="action: list or list")

args = parser.parse_args()

if args.filetype == 'genomes':
    map(print, genomes_list())

if args.filetype == 'blastdb':
    map(print, blastdb_list())

if args.filetype == 'protosp':
    map(print, protosp_list())

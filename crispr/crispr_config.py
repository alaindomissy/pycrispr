# import sys,os
# sys.path.insert(0,os.path.abspath(__file__+"/../.."))

from __future__ import absolute_import, division, print_function
# from __future__ import unicode_literals
import argparse
from .crispr.config import genomes_path, protosp_path, genomes_list, protosp_list

parser = argparse.ArgumentParser(description='blast a protospacers bed file to shortlist psoiible off-targets locations')

parser.add_argument('filetype',
                     help='filetype: genomes or protosp' )

parser.add_argument('action',
                     help="action: list or list")


args = parser.parse_args()

if args.filetype == 'genomes':
    map(print, genomes_list())

if args.filetype == 'protosp':
    map(print, protosp_list())

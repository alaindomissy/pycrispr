from __future__ import print_function
import argparse
from config import genomes_path, protosp_path, genomes_list, protosp_list

parser = argparse.ArgumentParser()

parser.add_argument('filetype',
                     help='filetype: genomes or protosp' )

parser.add_argument('action',
                     help="action: list or list")


args = parser.parse_args()

if args.filetype == 'genomes':
    map(print, genomes_list())

if args.filetype == 'protosp':
    map(print, protosp_list())

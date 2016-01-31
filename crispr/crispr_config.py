from __future__ import print_function
import argparse
from config import genomes_path, protosp_path, genomes_list, protosp_list

parser = argparse.ArgumentParser()

parser.add_argument('filetype',
                     help='filetype: prsp or scaf' )

parser.add_argument('action',
                     help="action: list or list")


args = parser.parse_args()

if args.filetype == 'scaf':
    map(print, genomes_list())

if args.filetype == 'prsp':
    map(print, protosp_list())

from __future__ import absolute_import, division, print_function    # , unicode_literals
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("stringtoecho", help="echo the string you use here")

parser.add_argument("square", help="display a square of a given number", type=int)

args = parser.parse_args()
print(args.stringtoecho)
print(args.square**2)

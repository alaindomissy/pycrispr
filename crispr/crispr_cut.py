from __future__ import print_function
try:
    from io import StringIO        # python3
except ImportError:
    from StringIO import StringIO  # python2
import argparse
from cut import cut_fastafile #,cut_unicodestring


parser = argparse.ArgumentParser()

parser.add_argument('fastafilepath',
                     help="a sequence file name to load, extension should be one of: fasta, fastq, fastq.gz or fastq.gz")

# parser.add_argument('unicodestring',
#                     help='a sequence as a unicode string',
#                     type=lambda s: unicode(s, 'utf8')
#                     )

args = parser.parse_args()


# print(cut_unicodestring(args.unicodestring))

# print(cut_fastafile(args.fastafilepath))
for bedline in cut_fastafile(args.fastafilepath):
    print(bedline, end='')

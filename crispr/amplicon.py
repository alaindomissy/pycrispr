########################################################################################################################
# AMPLICON
#
# API functions:  amplicon
#
########################################################################################################################
from __future__ import absolute_import, division, print_function   # , unicode_literals
import operator
import pickle
# from Bio.Seq import Seq
from Bio import SeqIO as seqio
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from .score import get_position, is_forward
from .stretch import run
from .config import ampliconlog



# AMPLICONS IO
##############

def save_amplicons(amplicons, filename, directory):
    """
    save amplicons
    :param guides:
    :param filename:
    :param directory:
    :return:
    """
    with  open(directory + filename + ".amplicons.pkl", "w") as output_pkl:
        pickle.dump(amplicons, output_pkl)
    return amplicons


def load_amplicons(filename, directory):
    """
    load amplicons
    :param filename:
    :param directory:
    :return:
    """
    with  open(directory + filename + ".amplicons.pkl", "r") as input_pkl:
        amplicons = pickle.load(input_pkl)
    return amplicons


class Amplicon:
    '''
    class with methods for:
    - designing pcr primers
    - describing the guides within it
    properties:

    .genome: must match BLAST database name for screening
    .substr_seqrec:

    .start: index in the list of all guides for first good guide inside amplicon
    .end: index in the list of all guides for first bad guide after amplicon
    .guides_count: number of (good) guides in amplicon  (=len(self.guides))

    .guides: list of (good) guides inside amplicon
    .fiveprimeabut:  bad guide along substr_seqrec that abut good guides on left side
    .threeprimeabut: bad guide along substr_seqrec that abut good guides on right side



    .left_outside: guide to the left of the outside
    .right_outside: guide to the right of the outside
    .left_inside:
    .right_inside:

    .max_length: permissible region distance
    .permissible_start: absolute numbering vs substr_seqrec
    .required_start: absolute numbering vs substr_seqrec

    .required_start_relative:

    .permissible_region
    '''

    def __init__(self, start, end, guide_list, genome, substr_seqrec, substr_offset):
        self.genome = genome
        self.substr_seqrec = substr_seqrec
        self.substr_offset = substr_offset
        self.start = start
        self.end = end
        self.guides_count = end - start
        self.guides = guide_list[start:end]
        ### TODO should this 2 be split up ?
        try:
            self.fiveprimeabut = guide_list[start - 1]
            self.threeprimeabut= guide_list[end]
        except:
            self.fiveprimeabut = "0"
            self.threeprimeabut = "None"

        # TODO what happens if no abut (except cases just above)
        self.max_length = get_position(self.threeprimeabut) - get_position(self.fiveprimeabut)

        # how far from desired amplicon, primers can be designed without bumping into the next (non-specific) primer
        self.left_outside = self.fiveprimeabut
        self.right_outside = self.threeprimeabut
        self.left_inside = self.guides[0]
        self.right_inside = self.guides[-1]
        self.permissible_start = get_position(self.left_outside) + ( 10 if is_forward(self.left_outside) else 1)
        self.required_start = get_position(self.left_inside) + (18 if is_forward(self.left_inside) else 14)
        self.permissible_end = get_position(self.right_outside) + ( 19 if is_forward(self.right_outside) else 10)
        self.required_end= get_position(self.right_inside) + (8 if is_forward(self.left_inside) else 2)

        # Bounds that need to be included in PCR product :
        self.required_start_relative = self.required_start - self.permissible_start
        self.required_end_relative = self.required_end - self.permissible_start


        self.permissible_region = substr_seqrec[self.permissible_start - substr_offset : self.permissible_end - substr_offset]
        # Set up some other convenience stuff:
        self.permissible_region.name = get_position(self.fiveprimeabut)
        self.permissible_region.id = get_position(self.fiveprimeabut)
        self.permissible_region.description = str(self.guides_count)
        self.permissible_region.seq.alphabet = IUPACAmbiguousDNA()

    def __str__(self):
        return \
        "Maximum permissible amplicon length: " + str(self.max_length) + "\n" + \
        "bp position of permissible amplicon start-end: " + str(self.permissible_start) + "-" + str(self.permissible_end) + "\n" + \
        "Number of guides contained: " + str(self.guides_count) + "\n" + \
        "Required priming subregion: " + str(self.required_start_relative) + "-" + str(self.required_end_relative) + "\n"


###################
# MAIN API FUNCTION
###################

def amplicon1(filename, directory, genome, threshold):
    """
    Takes a directory and filename
    from which to pickle load a list of potential guides complete with scores
    and finds "runs" of high-scoring guides which, if PCR-amplified as sub-regions,
    will generate a set of highly-specific sgRNAs.
    The Amplicon data structure stores the output information.
    :param filename:
    :param directory: working directory
    :param genome: an identifier for the genome assembly that the guides/scores
            are derived from (e.g. hg18, xl71)
    :param threshold:
    :return:a list of Amplicon instances, sorted by decreasing number of guides
    """

    guides_sorted, runs  = run(filename, directory, threshold)
    # print("len(guides_sorted)", len(guides_sorted))
    # print("guides_sorted", guides_sorted)

    ampliconlog('\nPRINT RUNS from amplicon')
    ampliconlog("============")
    ampliconlog("len(runs)", len(runs))
    ampliconlog("runs", runs)

    with open(directory + filename + '.fasta') as fasta_input:
        substr_seqrec = list(seqio.parse(fasta_input, 'fasta'))[0]

    substr_offset =  int(filename.split('_')[1].split('-')[0])

    amplicons = [Amplicon(start_end[0],  start_end[1], guides_sorted, genome, substr_seqrec, substr_offset) for start_end in runs]
    # print("len(amplicons)", len(amplicons))
    # print("amplicons", amplicons)

    # sorts by decreasing number of guides
    sorted_amplicons = sorted(amplicons, reverse=True, key=operator.attrgetter("guides_count"))
    # print("len(sorted_amplicons)", len(sorted_amplicons))
    # print("sorted_amplicons", sorted_amplicons)

    ampliconlog('\n*******************************************************************')
    ampliconlog('AMPLICONS')
    ampliconlog('*******************************************************************')
    for idx, ampl in enumerate(sorted_amplicons):
        ampliconlog("Amplicon", idx, ":")
        ampliconlog(ampl)

    ampliconlog('SAVE SORTED AMPLICONS')
    ampliconlog('---------------------')
    save_amplicons(sorted_amplicons, filename, directory)

    return sorted_amplicons

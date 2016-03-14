from __future__ import absolute_import, division, print_function
# from __future__ import unicode_literals
import copy
import operator
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA


def mask_sequence(item):
    '''
    Simply rewrites lowercase parts of a Seq object into Ns; Primer3 can be told
    to avoid Ns, but not lowercase sequence. (lowercase is often used in
    genomics data to indicate repetitive/non-complex sequence, so these
    regions are often poor choices for primer binding sites)
    '''
    masked_list = []
    for nucleotide in list(item.seq):
        if nucleotide.isupper():
            masked_list.append(nucleotide)
        else:
            masked_list.append("N")


    masked_string = "".join(masked_list)
    masked_sequence = copy.copy(item)
    masked_sequence.seq = Seq(masked_string, IUPACAmbiguousDNA())
    return masked_sequence


class Amplicon:
    '''
    A complete Amplicon for use in EATING; the class contains methods for designing
    amplifying primers and describing the guides within it.

    properties:

    .start, .end: count of runs along chromosome
    .fiveprimeabut, .threeprimeabut: count of guides along chromosome that abut good guides

    .guides: list of good guides in run
    .max_length: permissible region distance

    .guide_count: number of guides in list (=len(self.guides))

    .genome: genome; must match BLAST database name for screening
    .chromosome

    .left_outside: the id of the guide to the left of the outside
    .left_inside: ...

    .permissible_start: absolute numbering vs chromosome

    .required_start_absolute: absolute numbering vs chromosome
    .required_start_relative
    .permissible_region
    '''

    def __init__(self, start_end, guide_list, genome, chromosome):

        self.start = start_end[0]
        self.end = start_end[1]
        self.genome = genome
        self.chromosome = chromosome

        self.guide_count = self.end - self.start
        self.guides = guide_list[self.start:self.end]

        try:
            self.fiveprimeabut = guide_list[self.start - 1]
            self.threeprimeabut= guide_list[self.end]
        except:
            self.fiveprimeabut = "0"
            self.threeprimeabut = "None"

        self.max_length = int(self.threeprimeabut.name) - int(self.fiveprimeabut.name)


        # Figure out how far from the desired amplicon, primers can be designed
        # without bumping into the next (non-specific) primer
        self.left_outside = self.fiveprimeabut.id[-1]
        self.left_inside = self.guides[0].id[-1]

        if self.left_outside == "F" and self.left_inside == "R":
            self.permissible_start = int(self.fiveprimeabut.name) + 10
            self.required_start_absolute = int(self.guides[0].name) +14

        elif self.left_outside == "R" and self.left_inside == "R":
            self.permissible_start = int(self.fiveprimeabut.name) + 1
            self.required_start_absolute = int(self.guides[0].name) +14

        elif self.left_outside == "R" and self.left_inside == "F":
            self.permissible_start = int(self.fiveprimeabut.name) + 1
            self.required_start_absolute = int(self.guides[0].name) +18

        elif self.left_outside == "F" and self.left_inside == "F":
            self.permissible_start = int(self.fiveprimeabut.name) + 10
            self.required_start_absolute = int(self.guides[0].name) +18
        else:
            print("error on left")

        # (fiveprimeabuttingguide, threeprimeabuttingguide), end-start, scores_and_details[start:end]))
        #self.right_inside = item[2][-1][1].id[-1]
        self.right_inside = self.guides[-1].id[-1]
        self.right_outside = self.threeprimeabut.id[-1]

        if self.right_outside == "F" and self.right_inside == "R":
            self.permissible_end = int(self.threeprimeabut.name) + 19
            self.required_end_absolute = int(self.guides[-1].name) + 2

        elif self.right_outside == "R" and self.right_inside == "F":
            self.permissible_end = int(self.threeprimeabut.name) + 10
            self.required_end_absolute = int(self.guides[-1].name) + 8

        elif self.right_outside == "R" and self.right_inside == "R":
            self.permissible_end = int(self.threeprimeabut.name) + 10
            self.required_end_absolute = int(self.guides[-1].name) + 2

        elif self.right_outside == "F" and self.right_inside == "F":
            self.permissible_end = int(self.threeprimeabut.name) + 19
            self.required_end_absolute = int(self.guides[-1].name) + 8
        else:
            print("error on right")

        self.permissible_region = self.chromosome[self.permissible_start:self.permissible_end]

        # Bounds that need to be included in PCR product :
        self.required_start_relative = self.required_start_absolute-self.permissible_start
        self.required_end_relative = self.required_end_absolute - self.permissible_start

        # Set up some other convenience stuff:
        self.permissible_region.name =str(self.fiveprimeabut.name)
        self.permissible_region.id =str(self.fiveprimeabut.name)
        self.permissible_region.description=str(self.guide_count)
        self.permissible_region.seq.alphabet = IUPACAmbiguousDNA()

    def __str__(self):
        return \
        "Maximum permissible amplicon length: " + str(self.max_length) + "\n" + \
        "bp position of permissible amplicon start: " + str(self.permissible_start) + "\n" + \
        "Number of guides contained: " + str(self.guide_count) + "\n" + \
        "Required priming subregion: " + str(self.required_start_relative) + "-" + str(self.required_end_relative) + "\n"


def find_amplicon_start_end_s(guide_list, threshold):
    """

    :param guide_list: the list of guides in a region. This list must be:
            SeqRecord objects with :
            - their position within the region in their "name" attribute
            - their specificity score in their "annotations["score"]" attribute
            - their length obtained with the len function
    :param threshold: (default=94) the score cutoff below which a run will be "broken"
            Raising this leads to shorter runs (shorter PCR products) leading
            to fewer guides produced on average. Lowering it produces more PCR
            products, but of lower mean specificity.
    :return: a list of Amplicon instances, sorted by decreasing number of guides
    """

    # sort by increasing position within region
    guide_list = sorted(guide_list, reverse=False, key=lambda a: int(operator.attrgetter("name")(a)))

    run_starts=[]
    run_ends=[]
    prev_item_good = False
    for index, guide in enumerate(guide_list):
        this_item_good = guide.annotations["score"] >= threshold and len(guide) > 19
        if this_item_good and not prev_item_good :
                run_starts.append(index)
                prev_item_good = True
        if not this_item_good and prev_item_good:
            prev_item_good = False
            run_ends.append(index)
    run_start_end_s = [ (start, end) for (start, end) in zip(run_starts, run_ends) if end - start > 3]
    return run_start_end_s


def find_amplicons(guide_list, genome, chromosome, threshold=94):
    """
    Takes a list of potential guides complete with scores and finds "runs" of
    high-scoring guides which, if PCR-amplified as sub-regions, will generate a
    set of highly-specific sgRNAs. The Amplicon data structure stores the output
    information.
    :param guide_list:
    :param genome: an identifier for the genome assembly that the guides/scores
            are derived from (e.g. hg18, xl71)
    :param chromosome: an identifier for the chromosome or contig that the guides
            are derived from
    :param threshold:
    :return:
    """
    run_start_end_s = find_amplicon_start_end_s(guide_list, threshold)
    amplicons = [Amplicon(start_end, guide_list, genome, chromosome) for start_end in run_start_end_s]
    # sorts by decreasing number of guides
    sorted_amplicons = sorted(amplicons, reverse=True, key=operator.attrgetter("guide_count"))
    return sorted_amplicons

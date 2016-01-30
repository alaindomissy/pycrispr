
import copy
import operator
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA



# mask_sequence and screen_primer aren't directly called by the user.
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
    .length: permissible region distance

    .guidecount: number of guides in list (=len(self.guides))

    .genomename: genomename; must match BLAST database name for screening
    .chromosome

    .left_outside: the id of the guide to the left of the outside
    .left_inside: ...
    .permissible_start: absolute numbering vs chromosome
    .required_start_absolute: absolute numbering vs chromosome
    .required_start_relative
    .permissible_region
    '''
    def __init__(self, run, guidelist, genomename, chromosome):
        self.start = run[0]
        self.end = run[1]
        try:
            self.fiveprimeabut = guidelist[self.start-1]
            self.threeprimeabut= guidelist[self.end]
        except:
            self.fiveprimeabut = "0"
            self.threeprimeabut = "None"

        self.guides = guidelist[self.start:self.end]
        self.length = int(self.threeprimeabut.name)-int(self.fiveprimeabut.name)
        self.guidecount = self.end - self.start

        self.genomename = genomename
        self.chromosome = chromosome

        # Figure out how far from the desired amplicon
        #  primers can be designed without bumping into the next (non-specific) primer.
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
        self.permissible_region.description=str(self.guidecount)
        self.permissible_region.seq.alphabet = IUPACAmbiguousDNA()

    def __str__(self):
        return \
        "Maximum permissible amplicon length: " + str(self.length) + "\n" + \
        "bp position of permissible amplicon start: " + str(self.permissible_start) + "\n" + \
        "Number of guides contained: " + str(self.guidecount) + "\n" + \
        "Required priming subregion: " + str(self.required_start_relative) + "-" + str(self.required_end_relative) + "\n"





def find_amplicons(guidelist, genomename, chromosome, threshold=94):
    '''
    Takes a list of potential guides complete with scores and finds "runs" of
    high-scoring guides which, if PCR-amplified as sub-regions, will generate a
    set of highly-specific sgRNAs. The Amplicon data structure stores the output
    information.

    Arguments:
        guidelist: the list of guides in a region. This list must be:
            SeqRecord objects with their position within the region in their
            "name" attribute and their specificity score in their
            "annotations["score"]" attribute. This should be the format output
            by the EATING guide scoring functions.
        genomename: an identifier for the genome assembly that the guides/scores
            are derived from (e.g. hg18, xl71)
        chromosome: an identifier for the chromosome or contig that the guides
            are derived from
        threshold (default=94): the score cutoff below which a run will be "broken"
            Raising this leads to shorter runs (shorter PCR products) leading
            to fewer guides produced on average. Lowering it produces more PCR
            products, but of lower mean specificity.

        Returns a list of Amplicon instances.
    '''
    amplicon_positions = []
    guidelist = sorted(guidelist, reverse=False, key=lambda a: int(operator.attrgetter("name")(a)))
    scores = [(details.annotations["score"], len(details)) for details in guidelist]
    runs=[]
    ends=[]
    previtemgood = 0
    for index, (item, guide_length) in enumerate(scores):
        if item >= threshold and previtemgood ==0 and guide_length > 19:
            runs.append(index)
            previtemgood = 1
        elif item >= threshold and previtemgood == 1 and guide_length > 19:
            None
        elif previtemgood == 1:
            previtemgood = 0
            ends.append(index)
    runs = zip(runs, ends)

    for i in runs:
        if i[1] - i[0] > 3:
            amplicon_positions.append(Amplicon(i, guidelist, genomename, chromosome))
    goodruns = sorted(amplicon_positions, reverse=True, key=operator.attrgetter("guidecount"))# sorts by the number of guides in list
    return goodruns
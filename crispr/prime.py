########################################################################################################################
# PRIME
#
# API functions:  prime
#
########################################################################################################################
from __future__ import absolute_import, division, print_function    # , unicode_literals
import copy
import time
import primer3
from os.path import isfile
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from .config import PRIMER3_GLOBAL_ARGS
from .amplicon import amplicon
from .pcr import epcr_screen_primers, dumb_screen_primer
from .digest import nonoverlapping_guidecount
from .score import get_position
from .config import primelog


def mask_sequence(item):
    '''
    Rewrites lowercase parts of a Seq object into Ns;
    Primer3 can be told to avoid Ns, but not lowercase sequence.
    Lowercase is often used in genomics data to indicate repetitive/non-complex sequence,
    so these regions are often poor choices for primer binding sites.
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


def initialize_primers_list_file(filename):
    """
    Initialize a file with headers if it doesn't exist
    """
    headers_string ='\t'.join([
        "Sequence_id",
        "forward_seq", "forward_start", "forward_length", "forward_tm", "forward_gc",
        "reverse_seq", "reverse_start", "reverse_length", "reverse_tm", "reverse_gc",
        "input_seq_length",
        "PCR_product_length",
        "Guides_Contained",
        "Expanded priming distance",
        "Actual Non-overlapping Guide Count"])

    headers_string = headers_string + "\n"

    if not isfile(filename):
        with open(filename, "w") as handle:
            handle.write(headers_string)


def primer_search(ampl, filename, method="dumb", tm=40, primer3_global_args=PRIMER3_GLOBAL_ARGS):
    '''
    Returns a dict of candidate primers and their parameters, calculated against an Amplicon instance.

    In many cases, it's not possible to find perfect primers that cover all of the desired guides without
    also priming into some low-specificity guides.

    So, if this function can't find primers on the first round, it expands the permissible window at the
    possible expense of some guides within the amplicon.

    Depending on the primer pair chosen, the number of guides contained in a product will vary.
    The next function screens primers, and when the final pair is chosen the total number of primers is reported.

    Primer parameters can be supplied with the primer3_global_args argument, or left as defaults.
    '''
    #
    print("INITIALIZE AMPLIFICATIONPRINMERS.TSV FILE")
    initialize_primers_list_file(filename)

    # Work on a "masked" version of the sequence where lowercase letters are converted to Ns
    sequence_string = str(mask_sequence(ampl.permissible_region).seq)
    # print("sequence_string:" , sequence_string)

    increment_bases = 4 # number of nucleotides by which priming window is expanded at each round
    expand = 0 # Number of nucleotides current attempt shifts priming window by
    loopround = 0 # Round of primer attempts
    expansion_limit = ampl.max_length/8

    primerdict = {}
    primerdict["PRIMER_PAIR_NUM_RETURNED"] = 0
    i = 0
    is_bad = True
    genome = ampl.genome

    # Do this until a good primer is found.
    while is_bad and expand < expansion_limit:
        while primerdict["PRIMER_PAIR_NUM_RETURNED"] == 0 and expand < expansion_limit:
            expand = increment_bases * loopround
            # print("expand:" , expand)
            primeableregionleft_start = 0 + expand
            primeableregionleft_length = ampl.required_start_relative # + expand  # TODO remove the expand ?
            primeableregionright_start = ampl.required_end_relative - expand
            primeableregionright_length = len(ampl.permissible_region)-ampl.required_end_relative # + expand #TODO remove the expand ?
            primer3_seq_args = {
                "SEQUENCE_ID": ampl.guides[0].id,
                "SEQUENCE_TEMPLATE": sequence_string,
                "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": [[primeableregionleft_start, primeableregionleft_length,
                                                         primeableregionright_start,primeableregionright_length]],
                }
            primer3_global_args["PRIMER_PRODUCT_OPT_SIZE"] = int(ampl.max_length)
            # print("SEQUENCE_PRIMER_PAIR_OK_REGION_LIST:" , primer3_seq_args["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"])
            # Don't try priming if it's all Ns in the primeable region:
            # leftNs = sequence_string[primeableregionleft_start:primeableregionleft_start+primeableregionleft_length].
            # if

            ########################################################################
            primerdict = primer3.bindings.designPrimers(primer3_seq_args, primer3_global_args)
            print("\nPRIMER DICT:", primerdict)
            # prirmerdict = {}
            ########################################################################

            if primerdict["PRIMER_PAIR_NUM_RETURNED"] == 0:
                primelog(primerdict)
            primerdict["expandedpriming"] = expand
            loopround = loopround + 1
            primelog(str("Expanded Priming: " + str(expand) + " nt; " + str(primerdict["PRIMER_PAIR_NUM_RETURNED"]) + " primers tested."))

        # Loop through found primers and screen for specificity
        j = 0
        while j < primerdict["PRIMER_PAIR_NUM_RETURNED"] and is_bad:
            # primelog(j)
            leftprimer =  primerdict[str("PRIMER_LEFT_" + str(j) + "_SEQUENCE")]
            leftprimer_start =  str(primerdict[str("PRIMER_LEFT_"+ str(j))][0])
            leftprimer_length =  str(primerdict[str("PRIMER_LEFT_"+ str(j))][1])
            leftprimer_gc =  str(primerdict[str("PRIMER_LEFT_"+ str(j) + "_GC_PERCENT")])
            leftprimer_tm =  str(primerdict[str("PRIMER_LEFT_"+ str(j) + "_TM")])

            rightprimer = primerdict[str("PRIMER_RIGHT_" + str(j) + "_SEQUENCE")]
            rightprimer_start =  str(primerdict[str("PRIMER_RIGHT_"+ str(j))][0])
            rightprimer_length =  str(primerdict[str("PRIMER_RIGHT_"+ str(j))][1])
            rightprimer_gc =  str(primerdict[str("PRIMER_RIGHT_"+ str(j) + "_GC_PERCENT")])
            rightprimer_tm =  str(primerdict[str("PRIMER_RIGHT_"+ str(j) + "_TM")])
            primelog("Testing primer pair " + str(j) + " : " + leftprimer + " " + rightprimer)
            product_len = int(rightprimer_start) - int(leftprimer_start) #rightprimer_start is already right edge of primer

            if method == "dumb":
                is_bad = True
                ########################################################
                is_left_bad = dumb_screen_primer(leftprimer, genome)
                is_right_bad = dumb_screen_primer(rightprimer, genome)
                # is_left_bad = False
                # is_right_bad = False
                ##########################################################
                if not is_left_bad and not is_right_bad:
                    is_bad = False

                if is_left_bad:
                    # None
                    print("iteration" + str(i) + "left primer" + leftprimer + "is bad")
                if is_right_bad:
                    # None
                    print("iteration" + str(i) + "right primer" + rightprimer + "is bad")
                if is_bad and i == primerdict["PRIMER_PAIR_NUM_RETURNED"]:
                    with open(filename, "a") as handle:
                        #############################################################
                        handle.write("All the primers were bad for this amplicon!\n")
                        #############################################################

            if method == "epcr":
                #print(leftprimer, rightprimer)
                ###########################################################################################
                (is_bad, viableproducts) = epcr_screen_primers((leftprimer, rightprimer), genome, tm=tm)
                # print(viableproducts)
                ###########################################################################################

            j += 1

        if is_bad:
            primerdict["PRIMER_PAIR_NUM_RETURNED"] = 0
        else:
            # good primer pair(no off-target amplification) is found, add it to the output tsv
            print("\nampl.permissible_region", ampl.permissible_region)
            print("leftprimer_start", leftprimer_start)
            print("rightprimer_start", rightprimer_start)
            print("rightprimer_length", rightprimer_length)
            target = ampl.permissible_region.seq[int(leftprimer_start):int(rightprimer_start)+int(rightprimer_length)]

            product_nonoverlapping_guidecount = nonoverlapping_guidecount(target)
            tsv_list = [str(get_position(ampl.guides[0])),
                        leftprimer,leftprimer_start, leftprimer_length, leftprimer_tm, leftprimer_gc,
                        rightprimer, rightprimer_start, rightprimer_length, rightprimer_tm, rightprimer_gc,
                        str(ampl.max_length),
                        str(product_len),
                        str(ampl.guides_count),
                        str(primerdict["expandedpriming"]),
                        str(product_nonoverlapping_guidecount)]
            with open(filename, "a") as handle:
                #################################
                print("tsv_list", tsv_list)
                handle.write("\t".join(tsv_list)  + "\n")
                #################################

            ampl.primers = {"name": ampl.guides[0].name,
                            "leftprimer": leftprimer,
                            "leftprimer_start": leftprimer_start,
                            "leftprimer_length": leftprimer_length,
                            "leftprimer_tm": leftprimer_tm,
                            "leftprimer_gc": leftprimer_gc,
                            "rightprimer": rightprimer,
                            "rightprimer_start": rightprimer_start,
                            "rightprimer_length": rightprimer_length,
                            "rightprimer_tm": rightprimer_tm,
                            "rightprimer_gc": rightprimer_gc,
                            "template_length": ampl.max_length,
                            "product_length": product_len,
                            "guidecount": ampl.guides_count,
                            "expandedpriming": primerdict["expandedpriming"],
                            "product_nonoverlapping_guidecount": product_nonoverlapping_guidecount
                            }
            ########################################
            primelog("Primer added to output file.")
            ########################################

###################
# MAIN API FUNCTION
###################
# def prime(amplicon_list, ouputfilename = "datetime", method = "dumb", tm=40):

def prime(filename, directory, genome, threshold, method = "dumb", tm=40):
    """
    Design good (specific) primer pairs to pcr ampify a given list of amplicons
    :param ouputfilename:
    :param method:
    :param tm:
    :return:
    """
    amplicons_list = amplicon(filename, directory, genome, threshold)

    # if ouputfilename == "datetime":
    #     ouputfilename = str("primerlist_" + time.strftime("%Y%m%d-%H%M%S", time.localtime()) + ".tsv")
    ouputfilename = filename + '.amplificationprimers.csv'

    print('\nPRINT PRIMERS')
    # result = "\n".join([ "Amplicon %s :\n%s" % (index, ampl) for index, ampl in enumerate(amplicons_list)])
    # print(result)
    # primer_search(amplicons_list[0], ouputfilename, method, tm)

    for i, ampl in enumerate(amplicons_list):
        print("\nAmplicon %s :" % i)
        primer_search(ampl, ouputfilename, method, tm)

    # return result

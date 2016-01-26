from __future__  import print_function

import os
import copy
import operator
import time
# import re
# import sys
# import itertools
# from collections import Counter

# import primer3

import Bio
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline



from buffet.amplicon import Amplicon
from buffet.digest import nonoverlapping_guidecount


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
    probeyield = []
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


'''
Functions for primer design.
'''
# mask_sequence and screen_primer aren't directly called by the user.
def mask_sequence(item):
    '''
    Simply rewrites lowercase parts of a Seq object into Ns; Primer3 can be told
    to avoid Ns, but not lowercase sequence. (lowercase is often used in
    genomics data to indicate repetitive/non-complex sequence, so these
    regions are often poor choices for primer binding sites)
    '''
    s = []
    for nucleotide in list(item.seq):
        if nucleotide.isupper():
            s.append(nucleotide)
        else:
            s.append("N")
    s = "".join(s)
    masked_sequence = copy.copy(item)
    masked_sequence.seq = Seq(s, IUPACAmbiguousDNA())
    return masked_sequence


def screen_primer_dumb(primer, genomename):
    '''
    Input is a primer as a string.
    '''
    currfile = open("currprimer.fa", "w")
    currfile.write(">" + str(primer) + "\n")
    currfile.write(str(primer))
    currfile.close()
    blastn_cline = NcbiblastnCommandline(query="currprimer.fa", db=genomename, \
    task = "blastn-short",outfmt=5, out="primerblast.tmp", max_target_seqs=100, num_threads = 8, dust="no")
    blastn_cline
    result = blastn_cline()
    badprimer = 0
    os.remove("currprimer.fa")
    # Parse data
    result_handle = open("primerblast.tmp")
    blast_record = NCBIXML.read(result_handle) # if there were multiple queries, use NCBIXML.parse(result_handle)
    # How many matches are there with more than 14 or matching bases?
    match14 = 0
    for x in blast_record.alignments:
        for y in x.hsps:
            if y.positives > 14:
                match14 = match14 + 1
    match15 = 0
    for x in blast_record.alignments:
        for y in x.hsps:
            if y.positives > 15:
                match15 = match15 + 1
    # Set a cutoff of:
    if match14 > 40:
        badprimer = 1
    elif match15 > 30: #was 10
            badprimer = 1
    os.remove("primerblast.tmp")
    return badprimer




########
# Begin test set of functions to screen primers by proximity and orientation.

def blast_primer(primer_pair, genomename):
    p = []
    p.append(SeqRecord(seq = Seq(primer_pair[0], IUPACAmbiguousDNA()), id=str(str(primer_pair[0])), name = "F", description = ""))
    p.append(SeqRecord(seq = Seq(primer_pair[1], IUPACAmbiguousDNA()), id=str(str(primer_pair[1])), name = "R", description = ""))
    filename = "primer.fasta"
    SeqIO.write(p, filename, "fasta")
    blastn_cline = NcbiblastnCommandline(query=filename, db=genomename, \
    task = "blastn-short",outfmt=5, out=filename + ".blast", max_target_seqs=15, max_hsps=100, num_threads = 7, evalue = 10, dust="no")
    blastn_cline()
    result_handle = open(filename + ".blast")
    hits = NCBIXML.parse(result_handle)
    hits = [item for item in hits]
    return hits


def parse_primer_hits(hits, tm=40):
    priming_dict = {}
    for item in hits:
        for align in item.alignments:
            for spot in align.hsps:
                genomic_binding_site = Bio.Seq.reverse_complement(str(spot.sbjct))

                ###########################################################################################
                #calculated_tm = primer3.bindings.calcEndStability(str(item.query), genomic_binding_site).tm
                calculated_tm = 55
                ###########################################################################################

                if calculated_tm > tm: # Calculating 3' end stability
                    try:
                        priming_dict[align.title].append((spot.sbjct_start, spot.sbjct_end, item.query, genomic_binding_site, calculated_tm))
                    except:
                        priming_dict[align.title] = [(spot.sbjct_start, spot.sbjct_end, item.query, genomic_binding_site, calculated_tm)]
    return(priming_dict)


def screen_hits(priming_dict):
    '''
    Group in 30kb intervals.
        If more than 6 hits in a 30kb interval, trash primer pair
        If fewer than 6 hits, see if they are in antiparallel orientation and
        get product size
    '''
    viableproducts = []
    for chromosome_hits in priming_dict.iteritems():
        product_counter = 0
        #print chromosome_hits[0]
        #Get the intervals between binding sites on a scaffold
        if len(chromosome_hits[1]) > 1 and len(chromosome_hits[1])<=6:
            ### For each primer, go through and ask if any other primer is within 30000 nt of it
            # If so, make a tuple with the two partners' starts and orientations
            # If they're different, then add it to a viable product list
            for binding_sites in chromosome_hits[1]:
                query_start = binding_sites[0]
                query_end = binding_sites[1]
                for other_sites in chromosome_hits[1]:
                    subject_start = other_sites[0]
                    subject_end = other_sites[1]
                    if abs(query_start - subject_start) <= 30000  and abs(query_start - subject_start) > 0:
                        #print("\t\tuhoh")
                        query_orientation = query_start - query_end # First primer has same index as interval list position
                        subject_orientation = subject_start - subject_end # Second primer has index+1 of interval list position
                        # Test for the sign (=direction) of the start/end subtraction for primer, if it's the same for both there's no pcr product.
                        sameorientation = all(chromosome_hits >= 0 for chromosome_hits in (query_orientation, subject_orientation)) or all(chromosome_hits < 0 for chromosome_hits in (query_orientation, subject_orientation))
                        # print("\t" + chromosome_hits[0] + " " + str(query_start - subject_start))
                        # print("\t\t" + str(sameorientation))
                        if sameorientation == False:
                            product_counter += 1
            viableproducts.append((str(chromosome_hits[0]), product_counter))
        if len(chromosome_hits[1]) > 6:
            product_counter = 999
            viableproducts.append((str(chromosome_hits[0]), product_counter))
    print(viableproducts)
    return viableproducts


def screen_primer_in_silico_pcr(primer_pair, genomename, tm=40):
    '''
    Outputs a 1 if the primer is bad, a 0 if it's good. Based on
    primer proximity and orientation within the desired genome.
    '''
    bad = 0
    hits = blast_primer(primer_pair, genomename)
    priming_dict = parse_primer_hits(hits, tm=tm)
    viableproducts = screen_hits(priming_dict)
    #if len(viableproducts) > 1 or viableproducts[0][1] > 2: #if there's a hits on more than one scaffold or more than one hit on a scaffold
    if sum([item[1] for item in viableproducts]) > 2:
        bad = 1
    return (bad, viableproducts)

# Default parameters for primer3 primer design - can supply your own by copying
# this and modifying. These defaults were intended for NEB's Q5 HotStart
# polymerase.
global_parameters = {
        'PRIMER_NUM_RETURN': 20,
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_PICK_INTERNAL_OLIGO': 0,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 50.0, # was 56
        'PRIMER_MAX_TM': 75.0, # was 67
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MAX_POLY_X': 100,
        'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_WT_SELF_ANY': 1, # added
        'PRIMER_PAIR_WT_COMPL_ANY':1, #added
        'PRIMER_PAIR_WT_COMPL_END':1, #added
        'PRIMER_PAIR_WT_SELF_END':1, #added
        'PRIMER_PAIR_WT_SELF_ANY':1, #added
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
        'PRIMER_PRODUCT_SIZE_RANGE': [26, 18000],
        'PRIMER_MIN_THREE_PRIME_DISTANCE':3,
        'PRIMER_PAIR_MAX_DIFF_TM':5, #added
        'PRIMER_PAIR_WT_DIFF_TM':1, #added
        'PRIMER_PAIR_WT_PRODUCT_SIZE_LT': 0.5, #added
        }

# User functions:
def primer_search(current_amp, global_parameters=global_parameters, filename="primerlist.tsv", method="dumb", tm=40):
    '''
    Returns a dict of candidate primers and their parameters, calculated
    against an Amplicon instance. In many cases, it's not possible
    to find perfect primers that cover all of the desired guides without
    also priming into some low-specificity guides. So, if this function can't
    find primers on the first round, it expands the permissible window at the
    possible expense of some guides within the amplicon.

    Depending on the primer pair chosen, the number of guides contained in a
    product will vary. The next function screens primers, and when the final
    pair is chosen the total number of primers is reported.

    Primer parameters can be supplied with the global_parameters argument, or
    left as defaults.
    '''
    # Initialize a file with headers if it doesn't exist
    try:
        with open(filename) as file:
            None
    except IOError:
        with open(filename, "w") as primerlist:
            primerlist.write("Sequence_id\tforward_seq\tforward_start\tforward_length\tforward_tm\tforward_gc\treverse_seq\treverse_start\treverse_length\treverse_tm\treverse_gc\tinput_seq_length\tPCR_product_length\tGuides_Contained\tExpanded priming distance\tActual Non-overlapping Guide Count\n")
            primerlist.close()

    # Work on a "masked" version of the sequence where lowercase letters
    # are converted to Ns
    sequence_string = str(mask_sequence(current_amp.permissible_region).seq)
    increment_bases = 4 # number of nucleotides by which priming window is expanded at each round
    expand = 0 # Number of nucleotides current attempt shifts priming window by
    loopround = 0 # Round of primer attempts
    expansion_limit = current_amp.length/8
    primerdict = {}
    primerdict["PRIMER_PAIR_NUM_RETURNED"] = 0
    i = 0
    bad = 1
    genomename = current_amp.genomename
    # Do this until a good primer is found.
    while bad == 1 and expand < expansion_limit:
        while primerdict["PRIMER_PAIR_NUM_RETURNED"] == 0 and expand < expansion_limit:
            expand = increment_bases * loopround
            primeableregionleft_start = 0 + expand
            primeableregionleft_length = current_amp.required_start_relative #+ expand #remove?
            primeableregionright_start = current_amp.required_end_relative - expand
            primeableregionright_length = len(current_amp.permissible_region)-current_amp.required_end_relative #+ expand #remove?
            seq_args = {"SEQUENCE_ID": current_amp.guides[0].id,
                     "SEQUENCE_TEMPLATE": sequence_string,
                     "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": [[primeableregionleft_start,primeableregionleft_length,
                                                             primeableregionright_start,primeableregionright_length]],
                     }
            global_parameters["PRIMER_PRODUCT_OPT_SIZE"] = int(current_amp.length)
            #print(seq_args["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"])
            # Don't try priming if it's all Ns in the primeable region:
            # leftNs = sequence_string[primeableregionleft_start:primeableregionleft_start+primeableregionleft_length].
            # if


            ########################################################################
            # primerdict = primer3.bindings.designPrimers(seq_args, global_parameters)
            prirmerdict = {}
            ########################################################################

            if primerdict["PRIMER_PAIR_NUM_RETURNED"] == 0:
                print(primerdict)
            primerdict["expandedpriming"] = expand
            loopround = loopround + 1
            print(str("Expanded Priming: " + str(expand) + " nt; " + str(primerdict["PRIMER_PAIR_NUM_RETURNED"]) + " primers tested."))
        # Loop through found primers and screen for specificity
        j = 0
        while j < primerdict["PRIMER_PAIR_NUM_RETURNED"] and bad == 1:
            print(j)
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
            print("Testing " + leftprimer + " " + rightprimer)
            product_len = int(rightprimer_start) - int(leftprimer_start) #rightprimer_start is already right edge of primer
            if method == "dumb":
                bad = 1
                left_bad = screen_primer_dumb(leftprimer, genomename)
                right_bad = screen_primer_dumb(rightprimer, genomename)
                if left_bad == 0 and right_bad == 0:
                    bad = 0
                if left_bad ==1:
                    None
                    #print("iteration" + str(i) + "left primer" + leftprimer + "is bad")
                if right_bad == 1:
                    None
                    #print("iteration" + str(i) + "right primer" + rightprimer + "is bad")
                if bad == 1 and i == primerdict["PRIMER_PAIR_NUM_RETURNED"]:
                    with open(filename, "a") as primerlist:
                        primerlist.write("All the primers were bad for this amplicon!\n")
                        primerlist.close()

            if method == "in_silico_pcr":
                #print(leftprimer, rightprimer)
                (bad, viableproducts) = screen_primer_in_silico_pcr((leftprimer, rightprimer), genomename, tm=tm)
                #print(viableproducts)
            j = j + 1
        # If a good primer (no off-target amplification) is found,
        # add it to the output tsv
        if bad == 0:

            target = current_amp.permissible_region[int(leftprimer_start):int(rightprimer_start)+int(rightprimer_length)]
            product_nonoverlapping_guidecount = nonoverlapping_guidecount(target)

            with open(filename, "a") as primerlist:
                primerlist.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % \
                     (current_amp.guides[0].name,leftprimer,leftprimer_start,leftprimer_length,leftprimer_tm,leftprimer_gc,\
                     rightprimer,rightprimer_start,rightprimer_length,rightprimer_tm,rightprimer_gc,\
                     current_amp.length,str(product_len),current_amp.guidecount,str(primerdict["expandedpriming"]),str(product_nonoverlapping_guidecount)))
                primerlist.close()
            current_amp.primers = {"name": current_amp.guides[0].name,
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
                                "template_length": current_amp.length,
                                "product_length": product_len,
                                "guidecount": current_amp.guidecount,
                                "expandedpriming": primerdict["expandedpriming"],
                                "product_nonoverlapping_guidecount": product_nonoverlapping_guidecount
                                }
            print("Primer added to output file.")
        elif bad == 1:
            primerdict["PRIMER_PAIR_NUM_RETURNED"] = 0


def collect_good_primers(amplicon_list, filename = "datetime", method = "dumb", tm=40):
    if filename == "datetime":
        filename = str("primerlist_" + time.strftime("%Y%m%d-%H%M%S", time.localtime()) + ".tsv")
    for i, item in enumerate(amplicon_list):
        print("\nAmplicon " + str(i))
        primer_search(item, filename=filename, method=method, tm=tm)

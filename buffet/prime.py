from __future__  import print_function

import time

import primer3

from settings import PRIMER3_PARAMETERS
from amplicon import mask_sequence, Amplicon
from primers_screen import screen_primer_dumb, screen_primer_in_silico_pcr
from digest import nonoverlapping_guidecount


def primer_search(current_amp, global_parameters=PRIMER3_PARAMETERS, filename="primerlist.tsv", method="dumb", tm=40):
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

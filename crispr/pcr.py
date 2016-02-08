
from __future__ import absolute_import, division, print_function
# from __future__ import unicode_literals
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline


# TODO check emboss primersearch
# http://emboss.sourceforge.net/apps/release/6.2/emboss/apps/primersearch.html

# TODO check find_differential_primers
# https://github.com/widdowquinn/find_differential_primers

# TODO check BioPython Emboss/Primer3
# https://github.com/biopython/biopython/blob/master/Bio/Emboss/Applications.py

# TODO check BioStar question
# Question: Automated Primer-Blast like funcionnality
# https://www.biostars.org/p/107922/


# BLAST UTILITIES
#################
def blast_seqstring(seqstring, genome):
    with open("primer.fasta", "w") as handle:
        handle.write(">" + str(seqstring) + "\n")
        handle.write(str(seqstring))
    blastn_cline = NcbiblastnCommandline(query="primer.fasta",
                                         db=genome,
                                         task="blastn-short",
                                         outfmt=5,
                                         out="primer.blast",
                                         max_target_seqs=100,
                                         num_threads=4,
                                         dust="no")
    _ = blastn_cline()
    os.remove("primer.fasta")
    with open("primer.blast")as handle:
        blast_record = NCBIXML.read(handle)
    os.remove("primer.blast")
    return blast_record


def blast_primer_pair(primer_pair, genome):
    p = [SeqRecord(seq=Seq(primer_pair[0], IUPACAmbiguousDNA()), id=str(str(primer_pair[0])), name="F"),
         SeqRecord(seq=Seq(primer_pair[1], IUPACAmbiguousDNA()), id=str(str(primer_pair[1])), name="R")
         ]
    SeqIO.write(p, "primer.fasta", "fasta")
    blastn_cline = NcbiblastnCommandline(query="primer.fasta",
                                         db=genome,
                                         task="blastn-short",
                                         outfmt=5,
                                         out="prinmer.blast",
                                         max_target_seqs=15,
                                         max_hsps=100,
                                         num_threads=7,
                                         evalue=10,
                                         dust="no")
    _ = blastn_cline()
    # TODO can VCBIXML deal with a filemname string instead  of the handle?
    with open("primer.fasta") as handle:
        blast_records = NCBIXML.parse(handle)
    return list(blast_records)


# Sub functions for epcr_screen_pairs


def parse_primer_hits(hits, tm=40):
    priming_dict = {}
    for hit in hits:
        for alignment in hit.alignments:
            for hsp in alignment.hsps:
                genomic_binding_site = Seq.reverse_complement(str(hsp.sbjct))
                ###########################################################################################
                # CALL TO primer3
                #################
                # calculated_tm = primer3.bindings.calcEndStability(str(hit.query), genomic_binding_site).tm
                calculated_tm = 55
                is_3_prime_end_stable = calculated_tm > tm
                ###########################################################################################
                if is_3_prime_end_stable:
                    priming_dict.get(alignment.title, []).append(
                        (hsp.sbjct_start, hsp.sbjct_end, hit.query, genomic_binding_site, calculated_tm)
                    )
    return priming_dict


def screen_hits(priming_dict):
    """
    Screen primers by proximity and orientation
    Group in 30kb intervals.
        If more than 6 hits in a 30kb interval, trash primer pair
        If fewer than 6 hits, see if they are in antiparallel orientation and
        get product size

    :param priming_dict:
    :return:
    """
    viableproducts = []
    for chromosome_hits in priming_dict.iteritems():
        product_counter = 0
        # print chromosome_hits[0]
        # Get the intervals between binding sites on a scaffold
        if 1 < len(chromosome_hits[1]) <= 6:
                # For each primer, go through and ask if any other primer is within 30000 nt of it
                # If so, make a tuple with the two partners' starts and orientations
                # If they're different, then add it to a viable product list
                for binding_sites in chromosome_hits[1]:
                    query_start = binding_sites[0]
                    query_end = binding_sites[1]
                    for other_sites in chromosome_hits[1]:
                        subject_start = other_sites[0]
                        subject_end = other_sites[1]
                        if 0 < abs(query_start - subject_start) <= 30000:
                            # First primer has same index as interval list position
                            is_query_forward = query_start < query_end
                            # Second primer has index+1 of interval list position
                            is_subject_forward = subject_start - subject_end
                            same_orientation = is_query_forward == is_subject_forward
                            # if same orientation there's no pcr product.
                            if not same_orientation:
                                product_counter += 1
                viableproducts.append((str(chromosome_hits[0]), product_counter))
        if len(chromosome_hits[1]) > 6:
            product_counter = 999
            viableproducts.append((str(chromosome_hits[0]), product_counter))
    print(viableproducts)
    return viableproducts


####################
# OTHER API FUNCTION
####################
def dumb_screen_primer(primer_string, genome):
    blast_record = blast_seqstring(primer_string, genome)
    positives_values = [hsp.positives for alignment in blast_record.alignments for hsp in alignment.hsps]
    number_above_14 = sum([value > 14 for value in positives_values])
    number_above_15 = sum([value > 15  for value in positives_values])
    is_bad = number_above_14 > 40 or number_above_15 > 30
    return is_bad


###################
# MAIN API FUNCTION
###################
def epcr_screen_primers(primer_pair, genome, tm=40):
    """
    Outputs a 1 if the primer is bad, a 0 if it's good. Based on
    primer proximity and orientation within the desired genome.
    :param primer_pair:
    :param genome:
    :param tm:
    :return:
    """
    bad = 0
    hits = blast_primer_pair(primer_pair, genome)
    priming_dict = parse_primer_hits(hits, tm=tm)
    viableproducts = screen_hits(priming_dict)
    # if len(viableproducts) > 1 or viableproducts[0][1] > 2:
    # #if there's a hits on more than one scaffold or more than one hit on a scaffold
    if sum([product[1] for product in viableproducts]) > 2:
        bad = 1
    return bad, viableproducts

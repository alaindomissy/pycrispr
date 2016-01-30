import os
import Bio
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA


def screen_primer_dumb(primer, genomename):
    '''
    Input is a primer as a string.
    '''
    currfile = open("currprimer.fa", "w")
    currfile.write(">" + str(primer) + "\n")
    currfile.write(str(primer))
    currfile.close()
    blastn_cline = NcbiblastnCommandline(query="currprimer.fa",
                                         db=genomename,
                                         task = "blastn-short",
                                         outfmt=5,
                                         out="primerblast.tmp",
                                         max_target_seqs=100,
                                         num_threads = 8,
                                         dust="no")
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


# Screen primers by proximity and orientation
#############################################

def blast_primer(primer_pair, genomename):
    p = []
    p.append(SeqRecord(seq = Seq(primer_pair[0], IUPACAmbiguousDNA()), id=str(str(primer_pair[0])), name = "F", description = ""))
    p.append(SeqRecord(seq = Seq(primer_pair[1], IUPACAmbiguousDNA()), id=str(str(primer_pair[1])), name = "R", description = ""))
    filename = "primer.fasta"
    SeqIO.write(p, filename, "fasta")
    blastn_cline = NcbiblastnCommandline(query=filename,
                                         db=genomename,
                                         task = "blastn-short",
                                         outfmt=5,
                                         out=filename + ".blast",
                                         max_target_seqs=15,
                                         max_hsps=100,
                                         num_threads = 7,
                                         evalue = 10,
                                         dust="no")
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

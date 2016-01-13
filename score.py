########################################################################################################################
#
# SCORE
#
# API functions: score
#
########################################################################################################################


import subprocess
import cStringIO


from Bio.Blast import NCBIXML
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

# TODO this is also imported from cutbedlines.py, ok?
from Bio import SeqIO as seqio

from zhang import zhangscore


def is_valid_pam(pam):
    return True
    # return (len(pam)==3
    #     and (pam[1] == "G" or pam[1] == "A" or pam[1] == "g" or pam[1] == "a")
    #     and (pam[2] == "G" or pam[2] == "g"))



###################
# MAIN API FUNCTION
###################

def score(direct, fn_noext, blastdb_directory, blastdb_db, chunksize, nbrofchunks,
          reref_substrate_id=None, load_genome=False):

    # genomedict = a dict made from a FASTA file identical to the one used to make the BLAST DB.
    # dict keys should be BLAST db hit_def.
    if load_genome:
        print 'reading genome',
        genomeseq = seqio.parse(open(direct +'dict.fasta', 'rb'), "fasta", alphabet=IUPACAmbiguousDNA())
        genomedict = {}
        for seq in genomeseq:
            genomedict[seq.id] = seq
        print 'done'

    with open(direct + fn_noext + '.fasta') as guidesfn:
        guides = list(seqio.parse(guidesfn, "fasta"))

    blastrecords = []
    for chunknbr in range(1,nbrofchunks+1):
        fn_withext = fn_noext + str(chunksize)+ 'x' + str(chunknbr) + '.blast'
        blastfn = direct + fn_withext
        try:
            with open(blastfn) as blasthndl:
                print 'parsing', fn_withext,
                blastrecords.extend(list(NCBIXML.parse(blasthndl)))
                print 'done'
        except IOError as ioe:
            print 'skipping ', fn_withext

        # print blastrecords

    for blastindex, blastitem in enumerate(blastrecords):

        fullmatches = 0
        scorelist = []

        for alignment in blastitem.alignments:    # alignment corresponds to a hit in the blast file
                                                  # a hit is a whole seq from  blastdb, many hsps can exist for 1 hit
            for hsp in alignment.hsps:

                # getting the pam adjacent to the hsp's subject

                hit_threeprime_offset = len(guides[0]) - hsp.query_end
                pamstart = hsp.sbjct_end + hit_threeprime_offset

                # if sbjct_is_on_forward_strand, equivalent to (hsp.sbjct_end > hsp.sbjct_start)
                if (hsp.frame[1] > 0):
                    pam_start = hsp.sbjct_end + hit_threeprime_offset
                    pam_onebased_range = (pam_start + 1, pam_start + 3)
                    pam_zerobased_range = (pam_start , pam_start + 3)
                else:
                    pam_start = hsp.sbjct_end - hit_threeprime_offset
                    pam_onebased_range = (pam_start -3, pam_start -1)
                    pam_zerobased_range = (pam_start -4, pam_start -1)

                # TODO *** get the correct substrate id for subject (not same as hit!) - should be alignment.hit_def
                # TODO *** use betools and fai instead of loading full genome
                if load_genome:
                    # print 'using the genome dict to lookup pam'
                    lookup_context = genomedict[reref_substrate_id]
                    pam = lookup_context[pam_zerobased_range[0]:pam_zerobased_range[1]]
                else:
                    # print 'using blastdbcmd to lookup pam'
                    context_lookup_command = "blastdbcmd -db " + blastdb_directory + blastdb_db \
                                             + " -dbtype nucl -entry " + alignment.accession \
                                             + " -range %s-%s" % pam_onebased_range
                    context_lookup_process = subprocess.Popen(context_lookup_command, stdout=subprocess.PIPE, shell = True)
                    fstring = context_lookup_process.communicate()
                    fstring = cStringIO.StringIO(fstring[0])
                    pam = seqio.read(fstring, "fasta")

                if not (hsp.frame[1] > 0):
                    pam = pam.reverse_complement()

                # make match string padded to query length, where bar(|) is match and space( ) is non-match
                mmstr = list(hsp.match)
                if hsp.query_start > 1:
                    mmstr = list("." * (hsp.query_start - 1)) + mmstr
                if hsp.query_end < 20:
                    mmstr = mmstr + list("." * (20 - hsp.query_end))
                mmstr = "".join(mmstr)

                # TODO can this really not be the case? and then why do we not append a matchdit with score 0.0
                if len(pam) == 3:
                    matchdict = {"match_score": 0.0,
                                 "alignment": alignment.hit_def,
                                 "hit_location": hsp.sbjct_start,
                                 "hit_sequence": hsp.sbjct,
                                 "pam": str(pam.seq),
                                 "match_bars": mmstr}
                    if is_valid_pam(pam):  # Test for valid PAM adjacency
                        if (hsp.positives >16 and hsp.positives < 20):
                            matchdict["match_score"] = zhangscore(mmstr)

                        # "Zhang lab algorithm doesn't handle perfect matches: give it a 50 if it's perfect"
                        # I think we should give it a 100.0 !
                        # we dont want to count in scoring the first perfect match found,
                        # but we want to use all the others if any
                        if hsp.positives >= 20: #changed from == 20; might be worth keeping in mind for bugs
                            if fullmatches > 0:
                                matchdict["match_score"] = 100.0
                            fullmatches += 1
                    scorelist.append(matchdict)

        finalscore = int(10000.000 / (100.000 + float(sum(item["match_score"] for item in scorelist))))
        guides[blastindex].annotations['score'] = finalscore
        guides[blastindex].annotations["blastindex"] = blastindex

        # print "seq: ", guides[blastindex].seq, "score: ", finalscore, "id: ", guides[blastindex].id

        # print guides[blastindex].seq
        # print alignment.hit_def
        # print hsp.match
        # print hsp.sbjct
        # print pam.seq

    # print "HERE ARE THE GUIDES: "
    # print [(guide.id, guide.seq, guide.annotations.get('score')) for guide in guides]

    # sort guides by position, using the seqrecord id
    # guides.sort(key=lambda x: int(x.id.split(':')[2].split('-')[0]))
    guides.sort(key=lambda x: int(x.id.split(':')[1].split('-')[0]))

    return guides

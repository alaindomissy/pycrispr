########################################################################################################################
#
# SCORE
#
# API functions: score
#
########################################################################################################################


from __future__ import print_function
import subprocess
import cStringIO

from Bio.Blast import NCBIXML
from Bio.Application import ApplicationError
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio import SeqIO as seqio # TODO this is also imported from cut.py, ok?

# from buffet.settings import BLASTDBDIR
from buffet.zhang import zhangscore



def is_valid_pam(pam):
    return True
    # return (len(pam)==3
    #     and (pam[1] == "G" or pam[1] == "A" or pam[1] == "g" or pam[1] == "a")
    #     and (pam[2] == "G" or pam[2] == "g"))


def load_genome_dict(fastafilepath):
    # genomedict = a dict made from a FASTA file identical to the one used to make the BLAST DB.
    # dict keys should be BLAST db hit_def.
    genomeseq = seqio.parse(open(fastafilepath, 'rb'), "fasta", alphabet=IUPACAmbiguousDNA())
    genomedict = {}
    for seq in genomeseq:
        genomedict[seq.id] = seq
    return genomedict



###################
# MAIN API FUNCTION
###################

def score(direct, fn_noext, blastdb_db, chunk_size, nbrofchunks,
          reref_substrate_id=None, load_genome=False):

    if load_genome:
        fastafilepath = direct +'dict.fasta'
        print('\nLOADING REFERENCE GENOME from', end=' ')
        genomedict =load_genome_dict(fastafilepath)
        print('...done')
        print('will use genomedict for pam look up')
    else:
        print('will use blastdbcmdfor pam look up')

    print('start loading unscored guides', end='')
    with open(direct + fn_noext + '.fasta') as guidesfn:
        guides = list(seqio.parse(guidesfn, "fasta"))
        print(' ...done')

    print('\nPARSING BLAST RESULTS CHUNKS')
    blastrecords = []
    for chunknbr in range(1,nbrofchunks+1):
        fn_withext = fn_noext + '.' + str(chunk_size) + 'seqs.' + str(chunknbr) + '.blast'
        blastfn = direct + fn_withext
        try:
            with open(blastfn) as blasthndl:
                print('> parsing chunk ', str(chunknbr).zfill(3), fn_withext, end='')
                blastrecords.extend(list(NCBIXML.parse(blasthndl)))
                print(' done')
        except IOError as ioe:
            print('> missing chunk', chunknbr,  fn_withext, 'skipped')

        # print(blastrecords)

    for blastindex, blastitem in enumerate(blastrecords):

        fullmatches = 0
        scorelist = []
        print('>', 'blastitem:', blastitem)

        for alignment in blastitem.alignments:    # alignment corresponds to a hit in the blast file
                                                  # a hit is a whole seq from  blastdb, many hsps can exist for 1 hit
            print('  >', alignment.hsps[0].query, 'guide off-target-hits review')
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
                    print('  >', hsp.sbjct,'off-target-hit ref-genome pam-lookup' , end=' ')
                    print('+', end='')
                    lookup_context = genomedict[reref_substrate_id]
                    pam = lookup_context[pam_zerobased_range[0]:pam_zerobased_range[1]]
                else:
                    print('    >',  hsp.sbjct, 'off-target-hit blast-db pam-lookup', end=' ')
                    # print('*', end='')
                    fstring = ''
                    try:
                        context_lookup_command = "blastdbcmd -db " + blastdb_db \
                                             + " -dbtype nucl -entry " + alignment.accession \
                                             + " -range %s-%s" % pam_onebased_range
                        context_lookup_process = subprocess.Popen(context_lookup_command, stdout=subprocess.PIPE, shell = True)
                        fstring = context_lookup_process.communicate()
                        fstring = cStringIO.StringIO(fstring[0])
                    except ApplicationError as err:
                        print(str(err).split('message ')[1].strip('\''))
                    pam = seqio.read(fstring, "fasta") # if len(fstring)>0 else None

                if pam and not (hsp.frame[1] > 0):
                    pam = pam.reverse_complement()

                print('+ padding', end=' ')
                # make match string padded to query length, where bar(|) is match and space( ) is non-match
                mmstr = list(hsp.match)
                if hsp.query_start > 1:
                    mmstr = list("." * (hsp.query_start - 1)) + mmstr
                if hsp.query_end < 20:
                    mmstr = mmstr + list("." * (20 - hsp.query_end))
                mmstr = "".join(mmstr)

                print('+ pam-checking', end=' ')
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
                            print('+ zhang-scoring', end=' ')
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
                print()

        finalscore = int(10000.000 / (100.000 + float(sum(item["match_score"] for item in scorelist))))
        guides[blastindex].annotations['score'] = finalscore
        guides[blastindex].annotations["blastindex"] = blastindex
        print('> blastitem overall guide-scoring: ', finalscore, 'blastindex: ', blastindex)
        # print("seq: ", guides[blastindex].seq, "score: ", finalscore, "id: ", guides[blastindex].id)
        # print(guides[blastindex].seq)
        # print(alignment.hit_def)
        # print(hsp.match)
        # print(hsp.sbjct)
        # print(pam.seq)
    print(' done')


    print('\nCHECKING ZERO SCORED GUIDES')
    # TODO guides not getting a score, how does this happen ? fix it better
    # kind of fix guides without a score
    for guide in guides:
        try:
            guide.annotations['score']
        except:
            # print('guide %s does not have a score' % guide)
            guide.annotations['score'] = 0.0

    # print("HERE ARE THE GUIDES: ")
    # print([(guide.id, guide.seq, guide.annotations.get('score')) for guide in guides])

    print('\nSORTING SCORED GUIDES')
    # sort guides by position, using the seqrecord id
    # guides.sort(key=lambda x: int(x.id.split(':')[2].split('-')[0]))
    guides.sort(key=lambda x: int(x.id.split(':')[1].split('-')[0]))

    return guides

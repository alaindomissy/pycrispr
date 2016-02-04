########################################################################################################################
#
# SCORE
#
# API functions: score
#
########################################################################################################################


from __future__ import print_function
import subprocess
from Bio.Blast import NCBIXML
from Bio.Application import ApplicationError
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio import SeqIO as seqio # TODO this is also imported from cut.py, ok?
from config import scorelog
from zhang import zhangscore, format_factors
try:
    from io import StringIO        # python3
except ImportError:
    from StringIO import StringIO  # python2


def is_valid_pam(pam):
    return (len(pam)==3
        and (pam[1] == "G" or pam[1] == "A" or pam[1] == "g" or pam[1] == "a")
        and (pam[2] == "G" or pam[2] == "g"))


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

def score(nbrofchunks, filename, genome, direct, chunk_size, load_genome=False):

    filename = filename + '.prsp'
    blastdb = genome + '/' + genome

    if load_genome:
        fastafilepath = direct +'dict.fasta'
        print('\nLOAD REFERENCE GENOME from', end=' ')
        genomedict =load_genome_dict(fastafilepath)
        print('...done')
        scorelog('pam look up via genomedict')
    else:
        scorelog('pam look up via blastdbcmd')

    scorelog('load unscored guides', end='')
    with open(direct + filename + '.fasta') as guidesfn:
        guides = list(seqio.parse(guidesfn, "fasta"))
        scorelog(' ...done')

    print('\nPARSE BLAST-RESULTS IN CHUNKS')
    blastrecords = []
    for chunknbr in range(1, nbrofchunks + 1):
        fn_withext = filename + '.' + str(chunk_size) + 'seqs.' + str(chunknbr) + '.blast'
        blastfn = direct + fn_withext
        try:
            with open(blastfn) as blasthndl:
                print('\n> parse chunk ', str(chunknbr).zfill(3), fn_withext, end='')
                blastrecords.extend(list(NCBIXML.parse(blasthndl)))
                print(' done')
        except IOError as ioe:
            print('> miss chunk', chunknbr,  fn_withext, 'skip')

    # print(blastrecords)

    for blastindex, blastrecord in enumerate(blastrecords):

        fullmatches = 0
        scorelist = []
        # '(out of', nbrofchunks, ')   we dont know the number of prsps, although it is approx nbofchunks*chunk_size
        scorelog('  >', 'prsp', blastindex, 'hits', len(blastrecord.alignments), 'scaffold(s)')

        for alignment in blastrecord.alignments:  # alignment corresponds to a hit in the blast file
                                                  # a hit is a whole seq from  blastdb, many hsps can exist for 1 hit
            scorelog('    > ', alignment.title.split()[0], 'is hit with', len(alignment.hsps), 'hsps'),

            for hspnbr, hsp in enumerate(alignment.hsps):

                # getting the pam adjacent to the hsp's subject
                # hit_threeprime_offset = len(guides[0]) - hsp.query_end
                # TODO de-harcode the 20, above no longer working as guides i a generator now, not a list
                hit_threeprime_offset = 20 - hsp.query_end
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

                # TODO *** use betools and fai instead of loading full genome
                # if load_genome: DISABLED FOR NOW
                if False:
                    # scorelog('refgen', end=' ')
                    # TODO *** get the correct substrate id for subject (not same as hit!) - should be alignment.hit_def
                    substrate_id = 'TODO' #guide.id.split(':')[0]
                    lookup_context = genomedict[substrate_id]
                    pam = lookup_context[pam_zerobased_range[0]:pam_zerobased_range[1]]
                else:
                    # scorelog('blastdb', end=' ')
                    fstring = ''
                    try:
                        context_lookup_command = "blastdbcmd -db " + blastdb \
                                             + " -dbtype nucl -entry " + alignment.accession \
                                             + " -range %s-%s" % pam_onebased_range
                        context_lookup_process = subprocess.Popen(context_lookup_command, stdout=subprocess.PIPE, shell = True)
                        out, _ = context_lookup_process.communicate()
                        fstring = StringIO(out.decode())
                    except ApplicationError as err:
                        print(str(err).split('message ')[1].strip('\''))
                    pam = seqio.read(fstring, "fasta") # if len(fstring)>0 else None
                if pam and not (hsp.frame[1] > 0):
                    pam = pam.reverse_complement()
                # scorelog('pam', pam.seq, end=' ')


                # TODO can this really not be the case? and then why do we not append a matchdit with score 0.0
                # Test for valid PAM adjacency
                if len(pam) == 3:
                    if is_valid_pam(pam):

                        # scorelog('pam:', pam.seq, end=' ')
                        # scorelog('valid', end=' ')
                        # scorelog('+ padding', end=' ')
                        # make match string padded to query length, where bar(|) is match and space( ) is non-match
                        mmstr = list(hsp.match)
                        if hsp.query_start > 1:
                            mmstr = list("." * (hsp.query_start - 1)) + mmstr
                        if hsp.query_end < 20:
                            mmstr = mmstr + list("." * (20 - hsp.query_end))
                        mmstr = "".join(mmstr)
                        matchdict = {"match_score": 0.0,
                                     "match_factors": '',
                                     "alignment": alignment.hit_def,
                                     "hit_location": hsp.sbjct_start,
                                     "hit_sequence": hsp.sbjct,
                                     "pam": str(pam.seq),
                                     "match_bars": mmstr}
                        if hsp.positives <16:
                            matchdict["match_score"] = 0.0
                            matchdict["match_factors"] = '%s missmatches' % (20 - hsp.positives)
                        if hsp.positives >=16 and hsp.positives < 20:
                            ########################################################################
                            matchdict["match_score"], matchdict["match_factors"] = zhangscore(mmstr)
                            #########################################################################

                        # "Zhang lab algorithm doesn't handle perfect matches: give it a 50 if it's perfect"
                        # I think we should give it a 100.0 !
                        # we dont want to count in scoring the first perfect match found,
                        # but we want to use all the others if any
                        if hsp.positives >= 20: # TODO safe to chane tp ==20 ?
                            if fullmatches == 0:
                                matchdict["match_score"] = 0.0
                                matchdict["match_factors"] = 'first full match'
                            if fullmatches > 0:
                                matchdict["match_score"] = 100.0
                                matchdict["match_factors"] = format_factors(1, 1, 1, 100.0)
                            fullmatches += 1
                        scorelog('      %2s  sbjct %s ' % (hspnbr, hsp.sbjct), end ='')
                        scorelog('pam', pam.seq)
                        scorelog('          missm %s %s' % (mmstr, matchdict["match_factors"]))
                        scorelog('          query %s' % (hsp. query))

                        scorelist.append(matchdict)
                    else:
                        scorelog('      %2s  sbjct %s ' % (hspnbr, hsp.sbjct), end ='')
                        scorelog('nopam')
                        #scorelog('           squery ' % hspnbr, hsp.query, end ='')

        finalscore = int(10000.000 / (100.000 + float(sum(item["match_score"] for item in scorelist))))
        guides[blastindex].annotations['score'] = finalscore
        guides[blastindex].annotations['blastindex'] = blastindex

        guides[blastindex].description += ' blastindex %s prsp_score %s %%' % (blastindex, finalscore)
        # print(guides[blastindex].format('fasta'))

        # print(guides[blastindex],
        #       blastindex,guides[blastindex].
        #       annotations["blastindex"] ,
        #       guides[blastindex].annotations['score']
        #       )

        scorelog('  score: ', finalscore)
        # scorelog("seq: ", guides[blastindex].seq, "score: ", finalscore, "id: ", guides[blastindex].id)
        # scorelog(guides[blastindex].seq)
        # scorelog(alignment.hit_def)
        # scorelog(hsp.match)
        # scorelog(hsp.sbjct)
        # scorelog(pam.seq)


    print('\nCHECK ZERO-SCORED GUIDES')
    # TODO guides not getting a score, how does this happen ? fix it better
    # kind of fix guides without a score
    for guide in guides:
        try:
            guide.annotations['score']
        except:
            # scorelog('guide %s does not have a score' % guide)
            guide.annotations['score'] = 0.0

    # print("HERE ARE THE GUIDES: ")
    # print([(guide.id, guide.seq, guide.annotations.get('score')) for guide in guides])

    print('\nSORTING SCORED GUIDES')
    # sort guides by position, using the seqrecord id
    # guides.sort(key=lambda x: int(x.id.split(':')[2].split('-')[0]))

    guides.sort(key=lambda x: int( x.id.split(':')[1].split('-')[0]))

    with  open(direct + filename + ".scored.fasta", "w") as output_handle:
        seqio.write(guides, output_handle, "fasta")

    # print(guides.format('fasta'))

    # for guide in guides:
    #     print()
    #     # print(guide.format('gb'))   # TODO why this? Locus identifier 'chr:101153-101173(+)' is too long
    #     print(guide.format('fasta'))

    return guides

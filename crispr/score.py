########################################################################################################################
#
# SCORE
#
# API functions: score
#
########################################################################################################################

from __future__ import absolute_import, division, print_function
# from __future__ import unicode_literals
import subprocess
import pickle
from Bio.Blast import NCBIXML
from Bio.Application import ApplicationError
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio import SeqIO as seqio # TODO this is also imported from cut.py, ok?
from .config import scorelog, blastdb_path
from .cut import tabbed_string_from_list, savelinestofile
from .zhang import zhangscore, format_factors
try:
    from io import StringIO        # python3
except ImportError:
    from StringIO import StringIO  # python2


def is_valid_pam(pam):
    return (len(pam)==3
        and (pam[1] == "G" or pam[1] == "A" or pam[1] == "g" or pam[1] == "a")
        and (pam[2] == "G" or pam[2] == "g"))


def load_genome_dict(fastafilepath):
    """
    utility omly used by score function with option load_genome=True
    deprecated
    :param fastafilepath:
    :return: a dictionary made from a FASTA file identical to the one used to make the BLAST DB,
    the keys are BLAST db hit_def, the values are
    """
    genomeseq = seqio.parse(open(fastafilepath, 'rb'), "fasta", alphabet=IUPACAmbiguousDNA())
    genomedict = {}
    for seq in genomeseq:
        genomedict[seq.id] = seq
    return genomedict


# GUIDES UTILITIES
##################

"""
    >>>import Bio
    >>>guide1 = Bio.SeqRecord.SeqRecord(seq=Bio.Seq.Seq('ccatcccgtcaaagtcagcc', Bio.Alphabet.SingleLetterAlphabet()),
                    id='chr6:136640082-136640102(+)',
                    name='chr6:136640082-136640102(+)',
                    description='chr6:136640082-136640102(+) blastindex=20 score=100',
                    dbxrefs=[])
    >>>guide1.annotations["blastindex"]=20,
    >>>guide1.annotations["score"]=100,

    >>>is_forward(guide1)
    True
    >>>get_substrate(guide1)
    chr6
    >>>get_position(guide1)
    136640082
    >>>get_score(guide1)
    100 %
"""

def is_forward(guide):
    return guide.id.split('(')[1][0:1] == '+'

def get_substrate(guide):
    return guide.id.split(':')[0]

def get_score(guide):
    return guide.annotations['score']

def get_position(guide):
    return int(guide.name.split(':')[1].split('-')[0])

def get_score_s(guides):
    return [get_score(guide) for guide in guides]

def get_substrate_position_score_tuples(guides):
    return [(get_substrate(guide), get_position(guide), get_score(guide)) for guide in guides]


# GUIDES BED UTILITIES
######################

def scorebedline(guide, low, high):
    """
    output a guide as a one line bed format string, colored per low and high score thresholds
    :param guide:
    :param low:
    :param high:
    :return: a bed-format one line output string for given 'guide', with coloring based on 'low' and 'high' thresholds
    >>>import Bio
    >>>guide1 = SeqRecord(seq=Seq('ccatcccgtcaaagtcagcc', SingleLetterAlphabet()),
                    id='chr6:136640082-136640102(+)',
                    name='chr6:136640082-136640102(+)',
                    description='chr6:136640082-136640102(+) blastindex=20 score=100',
                    dbxrefs=[])
    >>>guide1.annotations["blastindex"]=20,
    >>>guide1.annotations["score"]=100,
    >>>scorebedline(guide1, 75, 75)

    """
    substrate_id = guide.id.split(':')[0]
    bedstart = guide.id.split(':')[1].split('-')[0]
    # rename scaffold field and make coord 'absolute' (ie based on higher level)
    # if reref_substrate_id:
    #     substrate_id = reref_substrate_id
    #     bedstart = str(int(bedstart) + int(guide.id.split(':')[1].split('-')[0]))
    bedend = str(int(bedstart) + 20)
    sense = guide.id[-2]
    try:
        score = str(guide.annotations['score'])
    except:
        score ='0'
    # TODO use gradual shades reflecting scores
    # suddenly got this error invalid literal for int() with base 10:0.0 Now changed int() to float() But why ?
    if float(score) < low:
        color = '255,0,0'       # red
    elif int(score) >= high:
        color = '0,255,0'       # green
    else:
        color = '0,0,255'       # blue
    cutinfo = [substrate_id, bedstart, bedend, 'crRNA', score, sense, bedstart, bedend, color]
    # cutinfo = [substrate_id, bedstart, bedend, 'crspr_guide', 'score', sense]
    bedline = tabbed_string_from_list(cutinfo)
    return bedline


def scorebedlines(guides, low, high):
    return [scorebedline(guide, low, high) for guide in guides]


# GUIDES IO
###########

def save_guides(guides, filename, directory, low, high):
    """
    save guides
    :param guides:
    :param filename:
    :param directory:
    :return:
    """
    with  open(directory + filename + ".guides.pkl", "w") as output_pkl:
        pickle.dump(guides, output_pkl)

    with  open(directory + filename + ".guides.fasta", "w") as output_fasta:
        seqio.write(guides, output_fasta, "fasta")

    savelinestofile(directory + filename + '.guides.bed',
                    scorebedlines(guides, low, high)
                    )

    return guides


def load_guides(filename, directory):
    """
    load guides
    :param filename:
    :param directory:
    :return:
    """
    with  open(directory + filename + ".guides.pkl", "r") as input_pkl:
        guides = pickle.load(input_pkl)
    # not necessary ?
    # sort by increasing position within region
    guides_sorted = sorted(guides, reverse=False, key=get_position)
    return guides_sorted


###################
# MAIN API FUNCTION
###################

def score(nbrofchunks, filename, genome, directory, chunk_size, load_genome=False, low=75, high=94,):
    """

    :param nbrofchunks:
    :param filename:
    :param genome:
    :param directory:
    :param chunk_size:
    :param load_genome:
    :return: guides a list of ...
    """

    # filename = filename + '.prsp'
    # blastdb = genome + '/' + genome
    blastdb = blastdb_path(genome)

    if load_genome:
        # TODO deprecate
        fastafilepath = directory +'dict.fasta'
        scorelog('\nLOAD REFERENCE GENOME from', end=' ')
        genomedict =load_genome_dict(fastafilepath)
        scorelog('...done')
        scorelog('pam look up via genomedict')
    else:
        scorelog('pam look up via blastdbcmd')

    scorelog('load unscored guides', end='')
    with open(directory + filename + '.prsp.fasta') as guidesfn:
        guides = list(seqio.parse(guidesfn, "fasta"))
        scorelog(' ...done')

    scorelog('\n*******************************************************************')
    scorelog('SCORE BLAST CHUNKS')
    scorelog('*******************************************************************')
    blastrecords = []
    for chunknbr in range(1, nbrofchunks + 1):
        fn_withext = filename + '.prsp.fasta.' + str(chunk_size) + 'seqs.' + str(chunknbr) + '.blast'
        blastfn = directory + fn_withext
        try:
            with open(blastfn) as blasthndl:
                scorelog('> parse chunk ', str(chunknbr).zfill(3), fn_withext, end='')
                blastrecords.extend(list(NCBIXML.parse(blasthndl)))
                scorelog(' done')
        except IOError as ioe:
            scorelog('> miss chunk', chunknbr,  fn_withext, 'skip')

    # scorelog(blastrecords)

    for blastindex, blastrecord in enumerate(blastrecords):

        fullmatches = 0
        scorelist = []
        # '(out of', nbrofchunks, ')   we dont know the number of prsps, although it is approx nbofchunks*chunk_size
        scorelog('  prsp%s' % blastindex, 'hits', len(blastrecord.alignments), 'scaffold(s)')

        for alignment in blastrecord.alignments:  # alignment corresponds to a hit in the blast file
                                                  # a hit is a whole seq from  blastdb, many hsps can exist for 1 hit
            scorelog('    ' + alignment.title.split()[0], 'with', len(alignment.hsps), 'hsps'),

            for hspnbr, hsp in enumerate(alignment.hsps):

                # getting the pam adjacent to the hsp's subject
                # hit_threeprime_offset = len(guides[0]) - hsp.query_end
                # TODO de-hardcode the 20, above no longer working as guides i a generator now, not a list
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
                        scorelog(str(err).split('message ')[1].strip('\''))
                    pam = seqio.read(fstring, "fasta") # if len(fstring)>0 else None
                if pam and not (hsp.frame[1] > 0):
                    pam = pam.reverse_complement()
                # scorelog('pam', pam.seq, end=' ')

                # make match string padded to query length, where bar(|) is match and space( ) is non-match
                # mmstr = list(hsp.match)
                leftpad = ''
                rightpad = ''
                if hsp.query_start > 1:
                    # leftpad = list("." * (hsp.query_start - 1))
                    leftpad = "." * (hsp.query_start - 1)
                if hsp.query_end < 20:
                    # rightpad =list("." * (20 - hsp.query_end))
                    rightpad = "." * (20 - hsp.query_end)
                mmstr = (leftpad + hsp.match + rightpad)

                # TODO can this really not be the case? and then why do we not append a matchdit with score 0.0
                # Test for valid PAM adjacency
                if len(pam) == 3:
                    if is_valid_pam(pam):
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
                        if hsp.positives >= 20: # TODO safe to chane to ==20 ?
                            if fullmatches == 0:
                                matchdict["match_score"] = 0.0
                                matchdict["match_factors"] = 'first full match'
                            if fullmatches > 0:
                                matchdict["match_score"] = 100.0
                                matchdict["match_factors"] = format_factors(1, 1, 1, 100.0)
                            fullmatches += 1
                        scorelog('      %2s  sbjct %s ' % (hspnbr, leftpad + hsp.sbjct + rightpad), end ='')
                        scorelog('pam', pam.seq)
                        scorelog('          missm %s %s' % (mmstr, matchdict["match_factors"]))
                        scorelog('          query %s' % (leftpad + hsp. query + rightpad))

                        scorelist.append(matchdict)
                    else:
                        scorelog('      %2s  sbjct %s ' % (hspnbr, leftpad + hsp.sbjct + rightpad), end ='')
                        scorelog('no pam')
                        #scorelog('           squery ' % hspnbr, hsp.query, end ='')

        finalscore = int(10000.000 / (100.000 + float(sum(item["match_score"] for item in scorelist))))
        guides[blastindex].annotations['score'] = finalscore
        guides[blastindex].annotations['blastindex'] = blastindex

        guides[blastindex].description += ' blastindex=%s score=%s' % (blastindex, finalscore)
        # scorelog(guides[blastindex].format('fasta'))

        # scorelog(guides[blastindex],
        #       blastindex,guides[blastindex].
        #       annotations["blastindex"] ,
        #       guides[blastindex].annotations['score']
        #       )

        scorelog('  prsp%s' % blastindex, 'score:', finalscore)
        # scorelog("seq: ", guides[blastindex].seq, "score: ", finalscore, "id: ", guides[blastindex].id)
        # scorelog(guides[blastindex].seq)
        # scorelog(alignment.hit_def)
        # scorelog(hsp.match)
        # scorelog(hsp.sbjct)
        # scorelog(pam.seq)


    scorelog('\nCHECK ZERO-SCORED GUIDES')
    # TODO guides not getting a score, how does this happen ? fix it better
    # kind of fix guides without a score
    for guide in guides:
        try:
            guide.annotations['score']
        except:
            # scorelog('guide %s does not have a score' % guide)
            guide.annotations['score'] = 0.0

    # scorelog("HERE ARE THE GUIDES: ")
    # scorelog([(guide.id, guide.seq, guide.annotations.get('score')) for guide in guides])

    scorelog('\nSORT SCORED GUIDES BY POSITION')
    # sort guides by position, using the seqrecord id
    # guides.sort(key=lambda x: int(x.id.split(':')[2].split('-')[0]))
    # guides.sort(key=get_position)

    scorelog('\nSAVE SORTED SCORED GUIDES')
    save_guides(guides, filename, directory, low, high)

    return guides

import protospacers

from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from Bio import Restriction as rst
amb = IUPACAmbiguousDNA()
seq1 = Seq('aaaaaccccctttttggggg', amb)
seq2 = seq1 + Seq(rst.HpaII.site, amb) + seq1 + Seq(rst.BfaI.site, amb) + seq1 + Seq(rst.ScrFI.site, amb) + seq1
seq3 = seq2 + seq2 + seq2
substrate3 = SeqRecord(seq3, id='substrate3')

print protospacers.create_analysis(seq3).full()
# {BfaI: [46, 139, 232], HpaII: [22, 115, 208], ScrFI: [71, 164, 257]})


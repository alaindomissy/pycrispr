import buffet.main as main

datadirpath = "/data/"

main.blast(datadirpath, 'chr6:136640001-136641000_1000.prsp', '/BLASTDB/', 'mm8')


from buffet import cut
DIRECT = '/data/scaffolds/mm8/'
cut.create_cutbedlines_from_seq_file(DIRECT, 'chr6-1-400', 'fasta')


from buffet import digest
DIRECT = '/data/scaffolds/hg38/'
digest.digest_file(DIRECT, 'hg38', 'fasta')


digest.digest_coord('./', 'chr21:42774475-42905100','mm8', 'fasta')


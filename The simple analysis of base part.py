import re

fasta = {}


with open('res\sequence.fasta') as file:
    for line in file:
        
        if line.startswith(">"):
            
            name = line[1:].rstrip()
            
            fasta[name] = ''
            continue
        
        fasta[name] += line.rstrip().upper()


print(fasta)

def nt_count(seq):
    ntCount = [0, 0, 0, 0]
    for gene in seq:
        for nt in gene:
            if nt == 'A':
                ntCount[0] += 1
            elif nt == 'T':
                ntCount[1] += 1
            elif nt == 'C':
                ntCount[2] += 1
            elif nt == 'G':
                ntCount[3] += 1
        nt_CG = format(((ntCount[2] + ntCount[3]) / len(gene)*100),'.6f')
    return ntCount, nt_CG

def re_RNA(seq):
    for seq_l in seq:
        RNAre = re.sub('T', 'U', seq_l)
    return RNAre

def rna_trans_protein(RNAre):
    codonTable = {
        'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
        'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
        'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
        'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
        'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
        'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
        'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
        'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
        'UAC':'Y', 'UAU':'Y', 'UAA':'', 'UAG':'',
        'UGC':'C', 'UGU':'C', 'UGA':'', 'UGG':'W',
    }
    proteinSeq = []

    proteinSeq = [codonTable[RNAre[i:i+3]] for i in range(0, len(RNAre), 3) if RNAre[i:i+3] in codonTable and codonTable[RNAre[i:i+3]] != '']
    return proteinSeq


Genes = list(fasta.values())

print(nt_count(Genes))

print(re_RNA(Genes))

print(rna_trans_protein(re_RNA(Genes)))
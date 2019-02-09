Gencode = {
    'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
    'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
    'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
    'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
    'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
    'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
    'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
    'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
    'UAC': 'Y', 'UAU': 'Y', 'UAA': 'STOP', 'UAG': 'STOP',
    'UGC': 'C', 'UGU': 'C', 'UGA': 'STOP', 'UGG': 'W'
}

def translate(file_name):
    with open(file_name, "r") as f:
        seq = f.read()
    amino_seq = ""
    flag = False
    i = 0
    while i <= len(seq) - 3:
        codon = seq[i: i + 3]
        if not flag and codon == "AUG":
            flag = True
        if flag:
            amino_seq += Gencode[codon]
            if codon in ('UAA', 'UGA', 'UAG'):
                break
            else:
                amino_seq += "-"
            i += 3
        else:
            i += 1
    return amino_seq


if __name__ == '__main__':
    prot_seq= translate("nuc1.txt")
    with open("protein.txt", "w") as f:
        f.write(prot_seq)


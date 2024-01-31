from collections import Counter

gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}



def is_valid_dna(seq : str) -> bool: # Checks if a DNA sequence contains only A, T, C, G nucleotides
    seq = seq.lower()
    for i in seq:
        if i not in ['a', 't', 'c', 'g']:
            return False
    return True

def get_gc_content(seq : str) -> int: # Returns the GC content of a DNA sequence
    if is_valid_dna(seq) == False:
        return -1
    seq = seq.lower()
    number_of_gc = seq.count("g") + seq.count("c")
    return number_of_gc / seq.__len__() * 100

def reverse_compliment(seq : str) -> str: # Returns the reverse compliment of a DNA sequence
    if is_valid_dna(seq) == False: # if seq is not a valid dna string returns empty string
        return ""
    seq = seq.upper()
    returned_seq = ""
    for i in seq:
        if i == "A":
            returned_seq = f'T{returned_seq}'
        if i == "T":
            returned_seq = f'A{returned_seq}'
        if i == 'C':
            returned_seq = f'G{returned_seq}'
        if i == 'G':
            returned_seq = f'C{returned_seq}'
    return returned_seq # return the reverse compliment of the sequence in caps.

def change_t_to_u(seq : str) -> str: # Changes T to U in a DNA sequence
    seq = seq.upper()
    seq = seq.replace("T","U")
    return seq

def get_transcription(seq : str, strand: int) -> str: # Returns the transcription of a DNA sequence
    if strand == -1:
        seq = reverse_compliment(seq)
    return change_t_to_u(seq)
    
def stats_amino_acids() -> dict:
    returned_dict : dict = {} 
    for value in gencode:
        if gencode[value] == '_':
            gencode[value] = "END"
        if gencode[value] in returned_dict:
            returned_dict[gencode[value]] = returned_dict[gencode[value]] + 1
        else:
            returned_dict[gencode[value]] = 1
            
    return returned_dict

def translate_seq(seq : str ,reading_frame : int = None ) -> list or str: # Returns the amino acid sequence of a DNA sequence
    seq = seq.upper()
    arrayed_seq = [seq[i:i+3] for i in range(0, len(seq), 3)]
    amino_acids = []
    print(arrayed_seq)
    if is_valid_dna(seq) == False:
        return ""
    if reading_frame == None:
        for i in range(3):
            # new_seq = seq[i:]
            # print(new_seq)
            amino_acids += translate_seq(seq, i + 1)

        return amino_acids
    else:
        new_seq = seq[reading_frame - 1:]
        arrayed_seq = [new_seq[i:i+3] for i in range(0, len(new_seq), 3)]
        for cur_seq in arrayed_seq:
            if len(cur_seq) == 3:
                amino_acids.append(gencode[cur_seq])
        return amino_acids

def count_codons(seq : str) -> dict: # Returns the count of codons in a DNA sequence
    seq = seq.upper()
    arrayed_seq = [seq[i:i+3] for i in range(0, len(seq), 3)]
    sumerized_dict = {}
    for i in arrayed_seq:
        if i in sumerized_dict:
            sumerized_dict[i] = sumerized_dict[i] + 1
        else:
            sumerized_dict[i] = 1
    return sumerized_dict


# print(gencode['AAA'])
print(stats_amino_acids()) # TODO send an email to check the stop by END if it is a "_" and that the question is not clear on what to do. 
seq = "ACGTTG"
seq = seq.lower()
print(get_gc_content(seq))
print(reverse_compliment(seq))
print(get_transcription(seq, -1))
print(translate_seq(seq))
print(count_codons(seq))
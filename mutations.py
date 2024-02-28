

 
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


baboon = "CTCCAGCCCCACCTGAACTCACCCCCAGGACCCCGCCACAAACCCCCCGTCTCCGCCTCCGGCGGGCCGAAGGTCAGGATGGGCCCTTCCACCTGCAATCTTCCCTTTTTCCCCTTCGCACGGTCAACCGCACGGTCAAACGCACGATCCAGTACTGGCCCTTTTCTGCCTCGGACCTCTATAATTGGAAAACCCATAACCCCTCCTTTTCCCAAGACCCCCAGGCCTTGACCTCGTTGATAGAATCAATTCTCCTCACCCACCAGCCTACCTGGGATGATTGTCAACAGCTTTTGCAGGTCCTTCTAACCACCGAAGAAAGGCAGCGAGTCCTCCTGGAAGCCCGGAAAAATCTGCCCCGGCCTGGAGGCCTTCCAACCCAGCTCCCCAATGAAATAGACGAGGGATTTCCCCTCACCCGCCCGGACTGGGATTATGAGACAGCACCGGGTAGGGAGAGTCTCCGAATCTATCGCCAGGCTCTGTTGGCGGGTCTCAAGGGGGCAGGAAAACGCCCCACCAATTTGGCCAAGGTAAGGACTATAACTCAGGGAAAAGATGAAAGCCCGGCAGCCTTTATGGAAAGACTTCTGGAAGGGTTTCGAATGTATACTCCATTTGATCCAGAAGCACCAGAACACAAGGCTACCGTGGCCATGTCGTTCATAGATCAGGCAGCACTAGACATAAAAGGAAAACTCCAAAGGCTAGACGGGATCCAAACTCATGGGCTGCAGGAATTAGTAAGGGAGGCAGAAAAGGTATACAATAAAAGGGAAACCCCAGAAGAAAGAGAAGCTAGGCTTATAAAAGAACAGGAAGAACGGGAAGTCGGAGAGACAGAAAAAGAGATAAGCATTTACCAAAATCTGGCAGCCGTAGTGACTGAAAAAAGGGCAGGAAAGTCAGGGGAAACAAGAAGGCGGCCTAAAGTAGATAAGGACCAGTGCGCCTACTGCAAAGAGCGAGGGCATTGGACCAAGGACTGCCCCAAGCCTCCTAGAGACCAGAAGAAACCCGCCCCTGTCCTCACCTTAGGTGAGGACAGCGAATAGGGGTGTCAGGGCTCTGGAGCCCCCCCCGAGCCCCGGCTAACTCTATCTGTAGGGGGGCATCCCACCACCTTCTTGGTGGACACAGGCGCCCAACACTCGGTTTTGACCAAGGCAAACGGACCCCTGTCCTCTCGTACATCTTGGGTCCAGGGGGCAACAGGAAGAAAGATGCACAAATGGACTAACCGCCGGACAGTTAACCTAGGGCAAGGAATGGTGACACACTCCTTCTTGGTGGTACCTGAATGTCCGTACCCCCTTCTGGGGCGAGATCTCCTAACCAAACTCGGAGCTCAGATCCACTTCTCCGAGGCAGGGGCCCAGGTGTTAGACCGAGATGGCCAACCCATCCAAATTTTGACTGTGTCTCTGCAAGATGAACACCGGCTTTTTGACATCCCGGTCACCACCAGCCTCC"
cat = "CTCCAGCCCCACCCGATCTAACTCCCCAGACCCCACCTCAGACCCCCCGCCTCCGTCTCCGGCGATCGGACGACCCGAACGGCCCTCCCACTTGGCAATCTTCCCTGTTTCCTCTCCGGACAGTCAACCGGACTGTCCAATATTGGCCATTTTCTGCATCCGATCTCTATAATTGGAAAACTCATAACCCTTCCTTTTCCCAAGAACCCCAGGCCCTAACCTCACTAATAGAGTCTATCCTCCTCACACACCAGCCTACCTGGGACGACTGCCAGCAGCTCTTACAGGTCCTTTTGACTACCGAGGAAAGGCAACGAGTTCTCCTGGAGGCCCGGAAAAATGTGCCGGGGCCTGGAGGTCTTCCAACTCAGCTTCCCAATGAAATAGATGAAGGGTTTCCCCTCACCCGCCCGGACTGGGACTATGAAACGGCTCCAGGTAGGGAGAGTCTCCGAATCTATCGCCAGGCTCTGTTGGCGGGTCTCAAGGGGGCAGGGAAACGGACCATAATCCAAGGAAAAGAGGAGAGCCCCGCAGCCTTTATGGAACGGCTTCTAGAGGGGTTTCGGATGTATACTCCATTTAACCCCGAGGCTCCAGAACATAAGGCTACCGTGGCCATGTCATTCATAGATCAGGCAGCGTCTGATATAAAAGGAAAGCTCCAACGGCTAGATGGGATTCAAACCTATGGGTTGCAGGAACTAGTTAGGGAGGCAGAAAAGGTATATAACAAAAGGGAGACCCCAGAGGAAAAGGAGGCTAGGCTAGCAAAAGAGCAGGAAGCACGGGAGGAACGAAGAGATCGAAAGAGAGATAAACATTTAACCAAAATCCTGGCAGCTGTAGTGACAGAAAACGGAACAGGAAAGTCAGGGGAAACGAAGAGGCGGCCCAAAGTAGGAAAAGACCAGTGCGCCTATTGCAAAGAACGCGGGCACTGGATCAAAGACTGCCCTAAGCGCCCCGGGAACTCAAAGAAACCTGCCTCTGTGCTCGCCTTAGGCGAAGAGAGCGAATAGGGATGTCAGGGCTCCGGAGCCCCCCCCGAGCCCCGGCTAACCCTATCTGTAGGGGGGCATCCCACCACCTTCCTGGTGGACACGGGAGCTCAACACTCGGTCTTAACCAAGGCAAATGGGCCTCTGTCCTCGCGTACCTCCTGGGTCCAAGGAGCGACAGGAAGAAAAATTCACAAATGGACTAACCGCCGCACGGTTGACTTAGGGCAAGGGACGGTGACACACTCCTTTCTGGTGGTACCCGAATGCCCGTACCCCCTTCTGGGACGAGACCTCCTAACCAAGCTTGGAGCCCAGATACACTTCTCCGAAACGGGGGCCCAAGTGCTAAACCGGGATGGCCAGCCTATCCAAGTCTTAACTGTGTCCTTACAAGACGAACACAGACTTTTTGAAACCCCGGTCACCACTAACCTCC"
macaque = "CCTCCCTCCACCCCATCTCCTTCTCCAGCCCCGCCCAAACTCACCCATCGGACGCCGCCACCGACACCTCCTCGTCTCCGCTTGCGGCGGACTGAGGACCCAGATGGCCCTTCCACTTGGCAATCCTCCCTTTTTCCCCTCTGTACCGTCAATCGCACGGTCCAGTACTGGCCCTTCTCTGCCTCTGACCTCTACAACTGGAAAACCCATAACCCTTCCTTTTCCCAAGACCCCCAGGCCCTAACCTCGTTGATAGAATCCATTCTCCTCACTCACCAGCCCACTTGGGATGATTGCCAGCAACTCTTGCAGGTCCTCCTAACCACTGAAGAAAGGCAGCGAGTCCTCCTGGAGGCCCGGAAAAATGTGCCAGGACCAGGAGGCCTCCCAACCCAACTTCCCAATGAAATAGACGAGGGATTTCCCCTCACCCGCCCGGACTAGGACTATGAAACGGCACCAGGTAGGGAGAGTCTCCGAATCTATCGCCAGGCTCTGTTGGCGGGTCTCAAAGGGGCAGGAAAGCACCCCACAAATTTGGCCAAGGTAAGGACCATAACTCAGGAAAGGGACGAAAGCCCAGCAGCCTTCATGGAAAGGCTTCTGCAAGGGTTCCGAATGTATACCCCATTCGATCCAGAGGCCCTAGAACATAAGGCTACCGTAGCTATGGCATTCATAGACCAAGCTGCATTAGATATCAAAGGAAAACTCCAAAGGCTAGATGGGATCCAAACCTATGGATTACAGGAATTGGTTAAGGAGGCAGAAAAAGTGTATAATAAGAGAGAAACCCCTGAGGGAAGGGAAGCCAGGTTAGCGAAGGAACAGATGGAGCGAGAGGATTGTAGGGACCGAGTGAGGGATAAGCATTTAACAAAAATCCTGGCGGCAGTTGTGAGAGAGAAAGGACCAGGGAGAGAGGGAGAGAAGCGGAGGCGGCCAAAAGTGGAAAAAGACCAGTGTGCCTACTGCAAAGAACGGGGACATTGGATCAAAGATTGCCCCAAGCGTCCTAAAGACCAAAAGAAACCTGCCGCTGTCCTCACCCTAGGTGAAGATAGTGAATAGGGGTGTCAAGGCTCTGGAGCCCCCCCCGAGCCCCAGCTAACTCTCTCTGTAGGGGGGCACCCCACCACCTTCTTGGTGGACACAGGGGCCCAACCTTCAGTTTTGACAAAGGCAGACGGGCCCGTGTCTTCCCGCACATCCTGGGTCCAGGGGGCAACAGGGGGAAAATTGCACAAGTGGACTAACCTCCGGACAGTTAATCTTGGACGAGGAATGGTGACACATTCCTTCTTGGTGGTACCTGAATGCCCATAGCCCCTCCTAGGGCGAGATCTTCTGACCAAGCTCGGAGCCCAAATCCATTTCTCCGAGACAGGGGCCCAGGTATTAAATTGGGACCGTCAGCCCATCCAAATCTTAACTGTGTCTCTGCAAGATGAATATCGGCTTTTTGACGCTCCAGTCATCACTAGCCTCC"


def check_for_mutations(correct_seq: str,seq: str,checked_mutation:str) -> bool:
    if seq == "":
        return False
    # seq = seq.upper()
    count = 0 
    min_len = min(len(seq),len(correct_seq))
    for j in range(min_len):
        i = 0
        # print(gencode[seq[j]])
        if gencode[seq[j]] == "_":
            return count
        while i < len(seq[j]):
            # print(correct_seq[j][i])
            if seq[j][i] == "C" and correct_seq[j][i] == "G" or seq[j][i] == "G" and correct_seq[j][i] == "C" and checked_mutation == "Nonsense":
                count += 1
            if correct_seq[j][i] == "A" and seq[j][i] == "T" or correct_seq[j][i] == "T" and seq[j][i] == "A" and checked_mutation == "Synonymous":
                count += 1
            i += 1
    return count
        
    
            
    # for i in seq:
    #     if i not in ['A', 'T', 'C', 'G']:
    #         count += 1
    return count


def is_valid_dna(seq : str) -> bool: # Checks if a DNA sequence contains only A, T, C, G nucleotides
    seq = seq.lower()
    for i in seq:
        if i not in ['a', 't', 'c', 'g']:
            return False
    return True

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



def translate_seq(seq : str ,reading_frame : int = None ) -> list or str: # Returns the amino acid sequence of a DNA sequence
    seq = seq.upper()
    arrayed_seq = [seq[i:i+3] for i in range(0, len(seq), 3)]
    amino_acids = []
    result_array = []
    # print(arrayed_seq)
    if is_valid_dna(seq) == False:
        return ""
    if reading_frame == None:
        # print(amino_acids)
        for i in range(3):
            # new_seq = seq[i:]
            # print(new_seq)
            new_amino = translate_seq(seq, i + 1)
            for j in new_amino:
                if len(j) == 3:
                    # print(j)
                    amino_acids += j
            # print(amino_acids)
        return amino_acids
    else:
        new_seq = seq[reading_frame - 1:]
        new_arrayed_seq = [new_seq[i:i+3] for i in range(0, len(new_seq), 3)]
        # print(new_arrayed_seq)
        for cur_seq in new_arrayed_seq:
            if len(cur_seq) == 3:
                # print(cur_seq)
                result_array.append(cur_seq)
        # print(result_array)
        return result_array

def count_mutation_by_type(position,type):
    baboon_sequenced = translate_seq(baboon,position)
    cat_sequenced = translate_seq(cat,position)
    maqueque_sequenced = translate_seq(macaque,position)
    if type == "nonsense":
        cat_baboon_differnce = check_for_mutations(baboon_sequenced,cat_sequenced,"Nonsense")
        baboon_macaque_differnce = check_for_mutations(baboon_sequenced,maqueque_sequenced,"Nonsense")
        print(f"The number of nonsense mutations in the cat is {cat_baboon_differnce}")
        print(f"The number of nonsense mutations in the macaque is {baboon_macaque_differnce}")
    if(type == "synonymous"):
        cat_baboon_differnce = check_for_mutations(baboon_sequenced,cat_sequenced,"Synonymous")
        baboon_macaque_differnce = check_for_mutations(baboon_sequenced,maqueque_sequenced,"Synonymous")
        print(f"The number of synonymous mutations in the cat is {cat_baboon_differnce}")
        print(f"The number of synonymous mutations in the macaque is {baboon_macaque_differnce}")
        

seq = "ACGTTG"
# seq = seq.lower()
# print(get_gc_content(seq))
# print(reverse_compliment(seq))
# print(get_transcription(seq, -1))
# print(translate_seq(baboon,1))
# baboon_translate = translate_seq(baboon,1)
# cat_translate = translate_seq(cat,1)
# macaque_translate = translate_seq(macaque,1)
count_mutation_by_type(3,"nonsense")
count_mutation_by_type(3,"synonymous")
print(translate_seq(seq,1))
# print(count_codons(seq))
def is_valid_dna(seq : str):
    seq = seq.lower()
    for i in seq:
        if i not in ['a', 't', 'c', 'g']:
            return False
    return True

def get_gc_content(seq : str):
    if is_valid_dna(seq) == False:
        return -1
    seq = seq.lower()
    number_of_gc = seq.count("g") + seq.count("c")
    return number_of_gc / seq.__len__() * 100

def reverse_compliment(seq : str):
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

def change_t_to_h(seq : str):
    seq = seq.upper()
    seq = seq.replace("T","H")
    return seq

def get_transcription(seq : str, strand: int):
    if strand == -1:
        seq = reverse_compliment(seq)
    return change_t_to_h(seq)
    
    
    
seq = "ACGTTG"
seq = seq.lower()
print(get_gc_content(seq))
print(reverse_compliment(seq))
print(get_transcription(seq, -1))
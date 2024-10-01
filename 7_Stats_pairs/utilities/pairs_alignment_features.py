from Bio.Align import substitution_matrices

def score_seqs(seq1, seq2, gap, substmat):
    """
    (Function adapted from Arnau Cordomí)

    >>> score_seqs("THEFASTCAT", "THEFATCAT ", 1, "BLOSUM62")
    29
    >>> score_seqs("THEFASTCAT", "THEFATCA T", 1, "BLOSUM62")
    34
    >>> score_seqs("THEFA TCAT", "THEFASTCAT", 1, "BLOSUM62")
    52
    >>> score_seqs("THEFASTCAT", "THE", 1, "BLOSSUM62")
    ''
    """
    #print(seq1, seq2)
    if seq1 == None or seq2 == None:
        return ""
    if len(seq1) == 0 or len(seq2) == 0:
        return ""
    if len(seq1) != len(seq2):
        return ""

    score = 0
    subst_mat = substitution_matrices.load(substmat)
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    for i in range(len(seq1)):
        if seq1[i] == "." or seq2[i] == "." or seq1[i] == '-' or seq2[i] == '-':
            score += gap

        else:
            try:
                score += subst_mat[seq1[i]][seq2[i]]
            except IndexError:
                continue

    return int(score)


def seq_identity(seq1, seq2):
    """
    Compute sequence identity between two sequences. We do not consider gaps present in same positions in both sequences.
    """
    if seq1 == None or seq2 == None:
        return ""
    if len(seq1) == 0 or len(seq2) == 0:
        return ""
    if len(seq1) != len(seq2):
        return ""

    equal_res = 0
    both_gap = 0
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    for i in range(len(seq1)):
        if seq1[i] == seq2[i] and seq1[i] != "." and seq1[i] != "-":
            equal_res += 1
        if seq1[i] == seq2[i] and (seq1[i] == "." or seq1[i] == "-"):
            both_gap += 1
        if (seq1[i] == "." and seq2[i] == "-") or (seq2[i] == "." and seq1[i] == "-"):
            both_gap += 1

    score = (100 * equal_res) / (len(seq1) - both_gap)
    score = round(score, 2)

    # FutureWarning: future pandas update will raise an error -> we show score in the float equivalent
    #score = f"{round(score, 2)}%"

    #print(seq1)
    #print('################# ################# #################')
    #print(seq2)
    #print(f'equal_res: {equal_res}')
    #print(f'both_gap: {both_gap}')
    #print(f'(100 * {equal_res}) / ({len(seq1)} - {both_gap})')
    #print(score)

    return score



def extract_sequences(pair, pfam_alignments):
    """
    Extract Pfam alignment sequences from both variants forming the pair.
    """
    uniprot1 = pair['uniprot1']
    # Save the number indicating the position of the aminoacidic change and transform it to float
    aa_change_pos1 = pair["pos_prot_change"]
    aa_change_pos1 = float(aa_change_pos1)

    uniprot2 = pair['uniprot2']
    
    #print(uniprot1, uniprot2)


    # Save the number indicating the position of the aminoacidic change and transform it to float
    aa_change_pos2 = pair["pos_prot_change2"]
    aa_change_pos2 = float(aa_change_pos2)
    
    if uniprot1:
        seq1 = search_correct_seq(pfam_alignments, uniprot1, aa_change_pos1)
    else:
        seq1 = ''
    if uniprot2:
        seq2 = search_correct_seq(pfam_alignments, uniprot2, aa_change_pos2)
    else:
        seq2 = ''
    return seq1, seq2


def search_correct_seq(file_content, uniprot, aa_change_pos):
    """
    Search correct sequence in Pfam file, corresponding to the desired UniProt code and the range of positions in which the variant relies.
    """
    # Iterate over the contents of the file
    for text in file_content:
        # Remove the newline characters
        text = text.replace("\n", "")
        if text.startswith(uniprot):
            # If yes, save the rang, initial and final position variables from the text
            rang = text.split("/")[1].split("-", maxsplit=1)
            pos_i = float(rang[0])
            pos_f = float(rang[1].split()[0])

            # Check if the aminoacidic change position fits inside the range of the alignment position
            if aa_change_pos >= pos_i and aa_change_pos <= pos_f:
                # Save the part of the text corresponding to the alignment
                alignment = text.split()[1]
                return alignment



def pairwise_alignment(seq1, seq2, return_alignment=False, gap_penalty=-2, mismatch_penalty = None, match_bonus = None):
    """
    Function return pairwise alignment or pairwise alignment score
    :param seq1: First seq
    :param seq2: Second seq
    :param return_alignment: Set this parameter to True if you want alignment, default false.
    :param gap_penalty: Gap penalty, default -2
    :return: Score or Alignment
    """
    # Initial initialization of first column and first row
    if not mismatch_penalty:
        mismatch_penalty = -1
    if not match_bonus:
        match_bonus = 1

    # Pointer_matrix
    left_pointer = 1
    top_pointer = 3
    top_left_pointer = 5

    seq1_len, seq2_len = len(seq1), len(seq2)

    score_matrix = [[0] * (seq1_len + 1) for _ in range(seq2_len + 1)]
    pointer_matrix = [[0] * (seq1_len + 1) for _ in range(seq2_len + 1)]
    for i in range(len(score_matrix)):
        score_matrix[i][0] = i * gap_penalty
        if i != 0:
            pointer_matrix[i][0] = top_pointer

    score_matrix[0] = [i * gap_penalty for i in range(seq1_len + 1)]
    for i in range(1, seq1_len + 1):
        pointer_matrix[0][i] = left_pointer

    for i in range(1, len(score_matrix)):
        for j in range(1, len(score_matrix[0])):
            top_score = score_matrix[i - 1][j] + gap_penalty
            left_score = score_matrix[i][j - 1] + gap_penalty
            if seq1[j - 1] != seq2[i - 1]:
                top_left_score = score_matrix[i - 1][j - 1] + mismatch_penalty
            else:
                top_left_score = score_matrix[i - 1][j - 1] + match_bonus
            max_score = max(top_score, left_score, top_left_score)
            score_matrix[i][j] = max_score

            if max_score == top_left_score:
                    pointer_matrix[i][j] += top_left_pointer
            elif max_score == top_score:
                    pointer_matrix[i][j] += top_pointer
            elif max_score == left_score:
                    pointer_matrix[i][j] += left_pointer

    final_alignment_score = score_matrix[len(score_matrix) - 1][len(score_matrix[0]) - 1]

    if return_alignment:
        alignment = get_alignment(seq1, seq2, pointer_matrix)
        return [alignment, final_alignment_score]

    return final_alignment_score


def get_alignment(seq1, seq2, arrow_matrix):
    """
    Function to return alignments based on the calculated scores
    :param seq1:
    :param seq2:
    :param arrow_matrix:
    :return:
    """
    ptr1, ptr2 = len(seq1) - 1, len(seq2) - 1
    final_alignment = [[], []]
    i = len(arrow_matrix) - 1
    j = len(arrow_matrix[0]) - 1
    while i > 0 and j > 0:
        if arrow_matrix[i][j] == 5:
            final_alignment[0].append(seq1[ptr1])
            final_alignment[1].append(seq2[ptr2])
            ptr1 -= 1
            ptr2 -= 1
            i -= 1
            j -= 1
        elif arrow_matrix[i][j] == 3:
            final_alignment[0].append("-")
            final_alignment[1].append(seq2[ptr2])
            ptr2 -= 1
            i -= 1
        elif arrow_matrix[i][j] == 1:
            final_alignment[0].append(seq1[ptr1])
            final_alignment[1].append("-")
            ptr1 -= 1
            j -= 1

    rev_alignment = [final_alignment[0][::-1], final_alignment[1][::-1]]
    return rev_alignment


a = "GCC"
b = "GTTCA"
print pairwise_alignment(a, b, return_alignment=True)

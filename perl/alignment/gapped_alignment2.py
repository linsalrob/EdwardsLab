__author__ = 'Rob Edwards'

import sys
from matrices import blosum62

def score(a, b):
    blosum = blosum62()
    if a in blosum and b in blosum[a]:
        return blosum[a][b]
    elif b in blosum and a in blosum[b]:
        return blosum[b][a]
    else:
        sys.stderr.write("Can not score amino acids " + a + " and " + b + "\n")
        return -8


def score_alignment(seq1, seq2, gap_open=11, gap_extend=1):
    """
    Generate a score for an alignment between two sequences. This does not do the alignment!
    :param seq1: The first sequence
    :param seq2: The second sequence
    :param gap_open: The gap opening penalty
    :param gap_extend: The gap extension penalty
    :return: An int for the best score for the alignment
    """

    score_matrix = [[[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)] for k in range(3)]

    # initially populate the gap open/extension
    for i in range(1, len(seq1)+1):
        score_matrix[0][i][0] = -gap_open - (i-1)*gap_extend
        score_matrix[1][i][0] = -gap_open - (i-1)*gap_extend
        score_matrix[2][i][0] = -10*gap_open

    for j in range(1, len(seq2)+1):
        score_matrix[2][0][j] = -gap_open - (j-1)*gap_extend
        score_matrix[1][0][j] = -gap_open - (j-1)*gap_extend
        score_matrix[0][0][j] = -10*gap_open

    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            lower_scores = [score_matrix[0][i-1][j] - gap_extend, score_matrix[1][i-1][j] - gap_open]
            score_matrix[0][i][j] = max(lower_scores)

            upper_scores = [score_matrix[2][i][j-1] - gap_extend, score_matrix[1][i][j-1] - gap_open]
            score_matrix[2][i][j] = max(upper_scores)

            mid_scores = [score_matrix[0][i][j], score_matrix[1][i-1][j-1] + score(seq1[i-1], seq2[j-1]), score_matrix[2][i][j]]
            score_matrix[1][i][j] = max(mid_scores)

    max_scores = [score_matrix[0][i][j], score_matrix[1][i][j], score_matrix[2][i][j]]
    return max(max_scores)


def gapped_alignment(seq1, seq2, gap_open=11, gap_extend=1):
    """
    Perform a gapped alignment. This approach uses two, 3 dimensional matrices.
    :param seq1: The first sequence
    :param seq2: The second sequence
    :param gap_open: The gap opening penalty (default = 11)
    :param gap_extn: The gap extention penalty (default = 1)
    :return: The score, and the two sequences with gaps in them
    """

    score_matrix = [[[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)] for k in range(3)]
    backtrack_matrix = [[[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)] for k in range(3)]

    # initially populate the gap open/extension
    for i in range(1, len(seq1)+1):
        score_matrix[0][i][0] = -gap_open - (i-1)*gap_extend
        score_matrix[1][i][0] = -gap_open - (i-1)*gap_extend
        score_matrix[2][i][0] = -10*gap_open

    for j in range(1, len(seq2)+1):
        score_matrix[2][0][j] = -gap_open - (j-1)*gap_extend
        score_matrix[1][0][j] = -gap_open - (j-1)*gap_extend
        score_matrix[0][0][j] = -10*gap_open

    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            lower_scores = [score_matrix[0][i-1][j] - gap_extend, score_matrix[1][i-1][j] - gap_open]
            score_matrix[0][i][j] = max(lower_scores)
            backtrack_matrix[0][i][j] = lower_scores.index(score_matrix[0][i][j])

            upper_scores = [score_matrix[2][i][j-1] - gap_extend, score_matrix[1][i][j-1] - gap_open]
            score_matrix[2][i][j] = max(upper_scores)
            backtrack_matrix[2][i][j] = upper_scores.index(score_matrix[2][i][j])

            mid_scores = [score_matrix[0][i][j], score_matrix[1][i-1][j-1] + score(seq1[i-1], seq2[j-1]), score_matrix[2][i][j]]
            score_matrix[1][i][j] = max(mid_scores)
            backtrack_matrix[1][i][j] = mid_scores.index(score_matrix[1][i][j])

    i=len(seq1)
    j=len(seq2)
    output_seq1 = seq1
    output_seq2 = seq2
    max_scores = [score_matrix[0][i][j], score_matrix[1][i][j], score_matrix[2][i][j]]
    max_score = max(max_scores)
    backtrack_level = max_scores.index(max_score)

    # we need this, time and again
    insert_indel = lambda word, i: word[:i] + '-' + word[i:]

    while i*j != 0:
        if backtrack_level == 0:
            if backtrack_matrix[0][i][j] == 1:
                backtrack_level = 1
            i -= 1
            output_seq2 = insert_indel(output_seq2, j)

        elif backtrack_level == 1:
            if backtrack_matrix[1][i][j] == 0:
                backtrack_level = 0
            elif backtrack_matrix[1][i][j] == 2:
                backtrack_level = 2
            else:
                i -= 1
                j -= 1
        else:
            if backtrack_matrix[2][i][j] == 1:
                backtrack_level = 1
            j -= 1
            output_seq1 = insert_indel(output_seq1, i)


    for k in xrange(i):
        output_seq2 = insert_indel(output_seq2, 0)
    for k in xrange(j):
        output_seq1 = insert_indel(output_seq1, 0)

    return (max_score, output_seq1, output_seq2)

if __name__ == "__main__":
    #s1 = 'ATGLVRRLGSFLVEDFSRYKLLL'
    #s2 = 'ATGLGLMRRSGSPLVESRYKLL'
    s1 = 'MQMCDRKHECYFEGFICDWHTLLEPHIVAQSEPYPCHKKMTQMPPPCSWFGNDIAEEKPSSIMATPAMPNVEEGM'
    s2 = 'MWMKDRKKNANECDWHPLLEYHIVAQSEPYKCCKKAMLGVKGAGTQMPPPCSWFGNDIAEEKPSSIMATPAMPNWEEGM'
    (score, s1, s2) = gapped_alignment(s1, s2)
    print(str(score) + "\n" + s1 + "\n" + s2)



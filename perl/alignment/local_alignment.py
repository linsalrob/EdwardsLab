__author__ = 'redwards'


import os
from array import array
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



def local_alignment(seq1, seq2, gap_open=11, gap_extn=1):
    """
    Perform a ungapped local alignment. This approach uses 6 matrices.
    :param seq1: The first sequence
    :param seq2: The second sequence
    :param gap_open: The gap opening penalty (default = 11)
    :param gap_extn: The gap extention penalty (default = 1)
    :return: The score, and the two sequences with gaps in them
    """

    ################################################
    # Create Score Matrix and Helper Matrices
    ################################################

    scorematrix = [array('i', [0 for x in xrange(0, len(seq1) + 1)]) for t in xrange(0, len(seq2) + 1)]
    maxhigheri = [array('i', [0 for x in xrange(0, len(seq1) + 1)]) for t in xrange(0, len(seq2) + 1)]
    maxhigherj = [array('i', [0 for x in xrange(0, len(seq1) + 1)]) for t in xrange(0, len(seq2) + 1)]
    helpermatrix = [array('i', [0 for x in xrange(0, len(seq1) + 1)]) for t in xrange(0, len(seq2) + 1)]

    endi = 0
    endj = 0
    bestscore = -1
    for i in xrange(0, len(seq2)):
        for j in xrange(0, len(seq1)):
            match = scorematrix[i][j] + score(seq2[i],seq1[j])
            maxhigheri[i+1][j+1] = max([maxhigheri[i+1][j] - gap_extn, scorematrix[i+1][j] - gap_open])
            maxhigherj[i+1][j+1] = max([maxhigherj[i][j+1] - gap_extn, scorematrix[i][j+1] - gap_open])
            maxv = max([match,maxhigheri[i+1][j+1],maxhigherj[i+1][j+1],0])
            scorematrix[i+1][j+1] = maxv

            # how did we get here?
            if maxv <= 0:
                helpermatrix[i+1][j+1] = 4
            elif maxv == match:
                helpermatrix[i+1][j+1] = 1
            elif maxv == maxhigherj[i+1][j+1]:
                helpermatrix[i+1][j+1] = 2
            elif maxv == maxhigheri[i+1][j+1]:
                helpermatrix[i+1][j+1] = 3
            else:
                helpermatrix[i+1][j+1] = 4

            newscore = scorematrix[i+1][j+1]
            if newscore > bestscore:
                bestscore = newscore
                endi = i+1
                endj = j+1



    i = endi
    j = endj
    hm = helpermatrix[i][j]
    while i * j != 0 and hm != 4:
        if hm == 1:
            i = i - 1
            j = j - 1
        elif hm == 2:
            i = i - 1
        elif hm == 3:
            j = j - 1

        hm = helpermatrix[i][j]

    print(bestscore)
    print(seq1[j:endj])
    print(seq2[i:endi])

if __name__ == "__main__":
    s1 = 'ATGLVRRLGSFLVEDFSRYKLL'
    s2 = 'GLMRRSGSPLVESRYKLL'
    local_alignment(s1, s2)

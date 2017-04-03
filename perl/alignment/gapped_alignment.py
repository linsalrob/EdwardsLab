__author__ = 'redwards'


import os
import sys
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



def gap_alignment(seq1, seq2, gap_open=11, gap_extn=1):
    """
    Perform a gapped alignment. This approach uses 6 matrices.
    :param seq1: The first sequence
    :param seq2: The second sequence
    :param gap_open: The gap opening penalty (default = 11)
    :param gap_extn: The gap extention penalty (default = 1)
    :return: The score, and the two sequences with gaps in them
    """

    scorematrix =   [array('i', [0 for x in xrange(0, len(seq1) + 1)]) for t in xrange(0, len(seq2) + 1)]
    insmatrix =     [array('i', [0 for x in xrange(0, len(seq1) + 1)]) for t in xrange(0, len(seq2) + 1)]
    inshelpmatrix = [array('i', [0 for x in xrange(0, len(seq1) + 1)]) for t in xrange(0, len(seq2) + 1)]
    delmatrix =     [array('i', [0 for x in xrange(0, len(seq1) + 1)]) for t in xrange(0, len(seq2) +1 )]
    delhelpmatrix = [array('i', [0 for x in xrange(0, len(seq1) + 1)]) for t in xrange(0, len(seq2) + 1)]
    helpermatrix =  [array('i', [0 for x in xrange(0, len(seq1) + 1)]) for t in xrange(0, len(seq2) + 1)]

    endi = 0
    endj = 0


    for i in xrange(1, len(seq2) + 1):
        scorematrix[i][0] = -gap_open - gap_extn*(i-1)
        insmatrix[i][0] = -gap_open - gap_extn*(i-1)
        delmatrix[i][0] = gap_open*(-10)#-d - de*(i-1)
    for j in xrange(1, len(seq1) + 1):
        scorematrix[0][j] = -gap_open - gap_extn*(j-1)
        insmatrix[0][j] = gap_open*(-10)#-d - de*(j-1)
        delmatrix[0][j] = -gap_open - gap_extn*(j-1)

    for i in xrange(0, len(seq2)):
        for j in xrange(0, len(seq1)):
            match = scorematrix[i][j] + score(seq2[i],seq1[j])

            inslist = [insmatrix[i][j+1] - gap_extn, scorematrix[i][j+1] - gap_open]
            insmatrix[i+1][j+1] = max(inslist)
            inshelpmatrix[i+1][j+1] = inslist.index(insmatrix[i+1][j+1])

            dellist = [delmatrix[i+1][j] - gap_extn, scorematrix[i+1][j] - gap_open]
            delmatrix[i+1][j+1] = max(dellist)
            delhelpmatrix[i+1][j+1] = dellist.index(delmatrix[i+1][j+1])

            scorelist = [match, insmatrix[i+1][j+1], delmatrix[i+1][j+1], 0]
            scorematrix[i+1][j+1] = max(scorelist)
            helpermatrix[i+1][j+1] = scorelist.index(scorematrix[i+1][j+1])



    i = len(seq2)
    j = len(seq1)

    news = ''
    newt = ''
    hm = 'helper'
    while i * j != 0:
        if hm == 'helper':
            if helpermatrix[i][j] == 1:
                hm = 'inshelp'
            elif helpermatrix[i][j] == 2:
                hm = 'delhelp'
            else:
                i = i - 1
                j = j - 1
                news = seq1[j] + news
                newt = seq2[i] + newt
        elif hm == 'inshelp':
            if inshelpmatrix[i][j] == 1:
                hm = 'helper'
            i = i - 1
            news = '-' + news
            newt = seq2[i] + newt
        elif hm == 'delhelp':
            if delhelpmatrix[i][j] == 1:
                hm = 'helper'
            j = j - 1
            news = seq1[j] + news
            newt = '-' + newt
        else:
            print 'should not get here'




    bestscore = scorematrix[len(seq2)][len(seq1)]
    return(bestscore, news, newt)


if __name__ == "__main__":
    # s1 = 'ATGLVRRLGSFLVEDFSRYKLLL'
    # s2 = 'ATGLGLMRRSGSPLVESRYKLL'
    s1 = 'MQMCDRKHECYFEGFICDWHTLLEPHIVAQSEPYPCHKKMTQMPPPCSWFGNDIAEEKPSSIMATPAMPNVEEGM'
    s2 = 'MWMKDRKKNANECDWHPLLEYHIVAQSEPYKCCKKAMLGVKGAGTQMPPPCSWFGNDIAEEKPSSIMATPAMPNWEEGM'
    gap_alignment(s1, s2)


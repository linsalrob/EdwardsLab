"""
Extract genotypes from metagenomes by concatenating matching sequences.
"""

import sys

import argparse
from roblib import sequences
import re

# data is going to be an array where each element is a hash
# the hash has to have the following elements:
#         seq :  this is the sequence and will be mutable - we will change the seqeunce
#         ids : this is a list of seqid that contribute to this sequence


data = []
longest = {'end':0, 'idx':-1} # which is the longest sequence we've seen so far so we can add to it in case of a gap
bases = {'A', 'T', 'G', 'C', 'N', '-'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Extract genotypes from metagenomes by concatenating matching sequences")
    parser.add_argument('-f', help='fasta sequence alignment file', required=True)
    parser.add_argument('-n', help='minimum number of sequences a genotype must be in (default = 1)', default=1, type=int)
    args = parser.parse_args()


    # to start we are just going to merge identical sequences
    byseq = {}
    for (seqid, seq) in sequences.stream_fasta(args.f):
        seq = seq.upper() # make sure we have only upper case
        seq = seq.strip() # strip of leading/trailing whitespace
        seq = seq.replace('N', '-')  # remove any N's in the sequence
        seq = seq.rstrip('-')  # remove all the trailing -s

        keep_seq = True
        for k in seq:
            if k not in bases:
                if keep_seq:
                    sys.stderr.write("Skipped {} becuase it has base {}\n".format(seqid, k))
                keep_seq = False
        if not keep_seq:
            continue

        if seq in byseq:
            byseq[seq].append(seqid)
        else:
            byseq[seq] = [seqid]

    # now we are going to sort the sequences by
    # start position (smallest to highest)
    # and then length (shortest to longest)
    # that way we should be adding a shorter sequence to a longer sequence

    allseqs = {}
    for seq in byseq:
        seqid = " ".join(byseq[seq])

        # calculate the start position and the length of the sequence
        m = re.match('(\-*)([ATGC\-]+)', seq)
        beg = len(m.group(1))
        seqlen = len(m.group(2))

        allseqs[seqid] = {'seq' : seq, 'startposition': beg, 'seqlen': seqlen}

    allids = sorted(allseqs, key=lambda x: allseqs[x]['startposition'] or allseqs[x]['seqlen'])

    for seqid in allids:
        seq = allseqs[seqid]['seq']

        if not data:
            data.append({'seq': seq, 'ids': [seqid]})
            longest['end'] = len(seq)
            longest['idx'] = 0
            continue

        # now find out where we should start
        beg = 0
        m = re.match('(\-*)([ATGC\-]+)', seq)
        beg = len(m.group(1))
        end = beg+len(m.group(2))

        # we want to see if we have this sequence as a consensus sequence from seq[beg:end] for any key in data
        added = False

        # check if there is a gap between the start of our sequence and end of all the other sequences
        if longest['end'] < beg:
            tempseq = data[longest['idx']]['seq']
            for i in (longest['end'], beg):
                tempseq += "-"
            tempseq += seq
            longest['end'] = end
            continue

        for i in range(len(data)):
            subseq = ""
            testseq = data[i]['seq']
            if len(testseq) < end:
                subseq = testseq[beg:]
            else:
                subseq = testseq[beg:end]

            # the order of these two is really important. we need to check first if our new
            # sequence is part of the old one. This also covers the case when they are the same!
            if m.group(2) in subseq:
                # our new sequence is a part of our old sequence. no need to added it
                break

            # check for a bidirectional partial match
            if subseq in m.group(2):
                # our new sequence is longer than our old sequence
                #sys.stderr.write("Adding 0-{}/{}-{}:\n{}\n{}\n".format(beg, beg, end, testseq, seq))
                if len(testseq) < beg:
                    testseq = testseq + seq[len(testseq):end]
                else:
                    testseq = testseq[0:beg] + seq[beg:end]
                sys.stderr.write("{}\n".format(testseq))
                data[i]['seq'] = testseq
                data[i]['ids'].append(seqid)
                added = True
                if end > longest['end']:
                    longest['end'] = end
                    longest['idx'] = i
                break



        if not added:
            sys.stderr.write("{} <<< {} <<< NEW\n".format(seq, seqid))
            data.append({'seq': seq, 'ids': [seqid]})
            if end > longest['end']:
                longest['end'] = end
                longest['idx'] = len(data) - 1


    seqcount = 0
    for i in data:
        if len(i['ids']) >= args.n:
            idtag = str(seqcount) +"_" + str(len(i['ids'])) + " [ids: " + " ".join(i['ids']) + "]"
            seqcount += 1
            print(">{}\n{}".format(idtag, i['seq']))
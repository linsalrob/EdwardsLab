import unittest

import sys

from rob import sequences
import os

class MyTestCase(unittest.TestCase):
    def test_read_fasta(self):
        # requires my dropbox files!
        fastafile = os.path.join(os.environ['HOME'], 'Dropbox/Metagenomics/51.hits_small.fa')
        gzipfile = fastafile + ".gz"
        lrzipfile = fastafile + ".lrz"

        if os.path.exists(fastafile):
            fa = sequences.read_fasta(fastafile)
            self.assertEqual(len(fa), 500)
            bps = [len(fa[i]) for i in fa]
            sumbp = sum(bps)
            self.assertEqual(sumbp, 54213)
            maxbp = max(bps)
            self.assertEqual(maxbp, 195)
        else:
            sys.stderr.write("WARNING: {} does not exist. Cannot test fasta files\n".format(fastafile))

        if os.path.exists(gzipfile):
            fa = sequences.read_fasta(gzipfile)
            self.assertEqual(len(fa), 500)
            bps = [len(fa[i]) for i in fa]
            sumbp = sum(bps)
            self.assertEqual(sumbp, 54213)
            maxbp = max(bps)
            self.assertEqual(maxbp, 195)
        else:
            sys.stderr.write("WARNING: {} does not exist. Cannot test fasta files\n".format(gzipfile))

        if os.path.exists(lrzipfile):
            fa = sequences.read_fasta(lrzipfile)
            self.assertEqual(len(fa), 500)
            bps = [len(fa[i]) for i in fa]
            sumbp = sum(bps)
            self.assertEqual(sumbp, 54213)
            maxbp = max(bps)
            self.assertEqual(maxbp, 195)
        else:
            sys.stderr.write("WARNING: {} does not exist. Cannot test fasta files\n".format(lrzipfile))
            return

if __name__ == '__main__':
    unittest.main()

"""
Merge multiple pdf files into
"""

import os
import sys
import argparse
import PyPDF2 as PDF

merger = PDF.PdfFileMerger(strict=False)
for f in sys.argv:
    if f.endswith('.pdf') and os.path.exists(f):
        sys.stderr.write("Adding {}\n".format(f))
        merger.append(PDF.PdfFileReader(f, 'rb'))
    else:
        sys.stderr.write("Skipped {}\n".format(f))

merger.write("AllDocs.pdf")

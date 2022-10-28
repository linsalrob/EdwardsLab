"""
Read the accession files (if they exist) and return a dict of accession.version -> taxid
"""
import sys
import os
import gzip


def read_acc_tax_id(dbtype: str, tax_dir: str, verbose: bool = False) -> dict[str, int]:
    """
    Read the accessions for either proteins or nucleotides and return a dict
    """

    data = {}
    if dbtype == 'prot':
        dbfile = 'prot.accession2taxid.FULL.gz'
        av_col = 0
        tid_col = 1
    elif dbtype == 'nucl':
        dbfile = 'nucl_wgs.accession2taxid.gz'
        av_col = 1
        tid_col = 2
    else:
        print(f"ERROR: Do not know what {dbtype} database type is!", file=sys.stderr)
        return data

    if os.path.exists(dbfile):
        dbfullfile = dbfile
    elif os.path.exists(os.path.join(tax_dir, dbfile)):
        dbfullfile = os.path.join(tax_dir, dbfile)
    elif os.path.exists(os.path.join(tax_dir, "accession2taxid", dbfile)):
        dbfullfile = os.path.join(tax_dir, "accession2taxid", dbfile)
    else:
        print(f"FATAL: We can not find {dbfile} anywhere!", file=sys.stderr)
        print(f"We looked in . {tax_dir}, and {tax_dir}/accession2taxid", file=sys.stderr)
        return data

    if verbose:
        print(f"Parsing the data in {dbfullfile}", file=sys.stderr)
    with gzip.open(dbfullfile, 'rt') as f:
        for l in f:
            p = l.strip().split("\t")
            data[p[av_col]] = p[tid_col]

    return data

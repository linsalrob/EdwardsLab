import os
import sys
from bs4 import BeautifulSoup
import argparse

import roblib

__author__ = 'Rob Edwards'


def parse_run(filename, verbose):
    """
    Parse a run.xml file. We mainly retrieve the SRA ID(s) from this file

    :param filename: the path and filename of the file to parse
    :type filename: str
    :param verbose: print more output
    :type verbose: bool
    :return: A set of SRA run IDs
    :rtype: set
    """

    sra_ids = set()
    if verbose:
        sys.stderr.write("Parsing run file: {}\n".format(filename))
    soup = BeautifulSoup(open(filename, 'r'), 'xml')
    for r in soup.RUN_SET.find_all("RUN"):
        for p in r.find_all("PRIMARY_ID"):
            sra_ids.add(p.text)

    return sra_ids


def parse_sample(filename, verbose):
    """
    Parse a sample.xml file. We mainly retrieve the metadata from this file

    :param filename: the path and filename of the file to parse
    :type filename: str
    :param verbose: print more output
    :type verbose: bool
    :return: A dict of metadata
    :rtype: dict
    :raises: KeyError
    """

    data = {}
    if verbose:
        sys.stderr.write("Parsing sample file: {}\n".format(filename))
    soup = BeautifulSoup(open(filename, 'r'), 'xml')
    for sample in soup.SAMPLE_SET.find_all("SAMPLE"):
        # get the identifier
        identifiers = sample.find("IDENTIFIERS")
        if identifiers:
            data['primary_id'] = identifiers.find("PRIMARY_ID").text
        else:
            raise KeyError("FATAL: no IDENTIFIERS tag found in {}".format(filename))

        # get the title
        title = sample.find("TITLE")
        if title:
            data['title'] = title.text

        # get the sample information
        si = sample.find('SAMPLE_NAME')
        if si:
            sin = si.find('SCIENTIFIC_NAME')
            if sin:
                data['scientific_name'] = sin.text
            sin = si.find('TAXON_ID')
            if sin:
                data['taxon_id'] = sin.text

        xrefs= []
        for sls in sample.find_all("SAMPLE_LINKS"):
            for sl in sls.find_all("SAMPLE_LINK"):
                for xr in sl.find_all("XREF_LINK"):
                    xrefs.append(xr.find("DB").text + "|" + xr.find("ID").text)
        data['xref'] = "; ".join(xrefs)

        for sas in sample.find_all("SAMPLE_ATTRIBUTES"):
            for sa in sas.find_all("SAMPLE_ATTRIBUTE"):
                tag = sa.find("TAG")
                val = sa.find("VALUE")
                if tag and val:
                    data[tag.text] = val.text
                elif tag:
                    sys.stderr.write("Found a tag {} for {} but no value\n".format(tag, filename))
                elif val:
                    sys.stderr.write("Found a value {} for {} but no tag\n".format(val, filename))

    return data


def parse_directory(directory, verbose):
    """
    Parse a directory and print out the data from the XML files therein

    :param directory: The path of the directory to parse
    :type directory: str
    :param verbose: print more output
    :type verbose: bool
    :return: dictionary of results
    :rtype: dict
    """

    if not os.path.exists(directory):
        raise IOError("FATAL: {} does not exist".format(directory))

    files = os.listdir(directory)
    runxmlfile = None
    samplefile = None
    # find the run.xml file
    for f in files:
        if f.endswith('run.xml'):
            if runxmlfile:
                sys.stderr.write("Crap, have two run.xml files\n")
            runxmlfile = f
        if f.endswith("sample.xml"):
            if samplefile:
                sys.stderr.write("Crap, have two sample.xml files\n")
            samplefile = f

    if not samplefile or not runxmlfile:
        return None

    data = {}
    if samplefile:
        data = parse_sample(os.path.join(directory, samplefile), verbose)
    else:
        sys.stderr.write("No sample.xml file for {}\n".format(directory))
    if runxmlfile:
        data['sra_ids'] = "; ".join(parse_run(os.path.join(directory, runxmlfile), verbose))
    else:
        sys.stderr.write("No run.xml file for {}\n".format(directory))
    return data


def parse_parent_directory(directory, verbose):
    """
    Parse the upper level parent directory of directories

    :param directory: directory name - the directory of directories
    :type directory: str
    :param verbose: print more output
    :type verbose: bool
    :return:
    :rtype:
    """

    if not os.path.exists(directory):
        raise IOError("FATAL: {} does not exist".format(directory))

    files = os.listdir(directory)

    results = {}
    tags = set()
    for f in files:
        if os.path.isdir(os.path.join(directory, f)):
            res = parse_directory(os.path.join(directory, f), verbose)
            if not res:
                continue
            results[f] = res
            for k in res:
                tags.add(k)
        else:
            sys.stderr.write("Skipped {} because it is not a directory\n".format(f))

    alltags = ['primary_id', 'title', 'scientific_name', 'taxon_id', 'xref', 'sra_ids']
    tags.difference_update(set(alltags))
    [alltags.append(x) for x in sorted(tags)]
    alltags = [roblib.ascii_clean(x) for x in alltags]
    print("DIRECTORY\t"+"\t".join(alltags))
    for readid in results:
        sys.stdout.write(readid)
        for t in alltags:
            sys.stdout.write("\t" + str(roblib.ascii_clean(results[readid].get(t, ""))))
        sys.stdout.write("\n")



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse a directory of directories of metadata. You can download the metadata tarball from ftp://ftp-trace.ncbi.nlm.nih.gov/sra/reports/Metadata/')
    parser.add_argument('-d', help='directory of directories of XML files', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    parse_parent_directory(args.d, args.v)

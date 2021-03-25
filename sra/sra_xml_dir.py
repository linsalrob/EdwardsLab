"""
Parse the NCBI XML and extract useful information!

In this version we
1. Generate a list of attributes as we go along
2. Only print out those attributes for which we have >1  entry. There are lots of attributes that appear to be unique to a single record, so we ignore those

"""

import os
import sys
import argparse
import xml.etree.ElementTree as ET
from roblib import message

# we make a list so that the order is guaranteed, and then make a set for O(1) lookup
# but later we just resort the sets :)
known_attrs = ['id', 'accession', 'last_update', 'access', 'publication_date', 'submission_date']
known_titles = ['Ids', 'Description - Title', 'Description - Comment', 'Owner - Name', 'Owner - Email', 'Release Date',
                'Links']

known_attrs_set = set(known_attrs)
known_titles_set = set(known_titles)
known_attributes_set = set()


def parse_ids(this_sample_id, ids, verbose=False):
    """
    Parse the IDs field
    :param this_sample_id: The current sample ID for debugging
    :param ids: The ID xml
    :param verbose: more output
    :return:
    """
    id_text = ""
    for anid in ids:
        if 'db' in anid.attrib:
            #dbs[anid.attrib['db']] = anid.text
            id_text = id_text + "{}|{}; ".format(anid.attrib['db'], anid.text)
        elif 'db_label' in anid.attrib and "Sample name" == anid.attrib['db_label']:
            id_text = id_text + "Sample|{}; ".format(anid.text)
        else:
            sys.stderr.write("WARNING: ({}) Id: {} text: {} is not a DB\n".format(this_sample_id, anid.tag, anid.text))
    if verbose:
        message(f"ID is {id_text}", 'blue')

    return id_text


def parse_links(this_sample_id, links, verbose=False):
    """
    Parse the links field
    :param this_sample_id:
    :param links:
    :param verbose: more output
    :return:
    """

    link_text=""
    for child in links:
        if 'type' in child.attrib and "entrez" == child.attrib['type']:
            link_text = link_text + "{}|{}; ".format(child.attrib['target'], child.text)
        elif 'label' in child.attrib and 'type' in child.attrib and "url" == child.attrib['type']:
            link_text = link_text + "url|[{}]({}); ".format(child.text, child.attrib['label'])
        elif 'type' in child.attrib and "url" == child.attrib['type']:
            link_text = link_text + "url|[{}]({}); ".format(child.text, child.text)

    return link_text

def parse_description(this_sample_id, desc, verbose=False):
    """
    Parse the description field and return the title and the
    :param this_sample_id: The current sample ID for debugging
    :param desc: The description xml
    :param verbose: more output
    :return:
    """
    desc_title = ""
    desc_comment = ""
    desc_org = ""

    for child in desc:
        if child.tag == 'Title':
            desc_title = child.text
        elif child.tag == 'Comment':
            desc_comment = child.text
        elif child.tag == 'Organism':
            if 'taxonomy_name' in child.attrib:
                desc_org = child.attrib['taxonomy_name']
        else:
            sys.stderr.write("In {} description we had a tag ({}) that we didn't parse\n".format(this_sample_id, child.tag))

    if desc_title:
        desc_title = desc_title.strip()

    if desc_comment:
        desc_comment = desc_comment.strip()

    if desc_org:
        desc_org = desc_org.strip()

    if verbose:
        message(f"Parsed description {desc_title}", "BLUE")

    return (desc_title, desc_comment, desc_org)


def parse_attribute(this_sample_id, attr, verbose=False):
    """
    Parse out the attribute field
    :param this_sample_id:
    :param attr:
    :param verbose: more output
    :return:
    """

    if 'harmonized_name' in attr.attrib:
        return attr.attrib['harmonized_name'], attr.text
    elif 'attribute_name' in attr.attrib:
        return attr.attrib['attribute_name'].lower(), attr.text
    else:
        sys.stderr.write("WARNING: ({}): Neither harmonized_name nor attribute name in attribute {}".format(this_sample_id, attr.text))
        return "",""


def parse_owner(this_sample_id, owner, verbose=False):
    """
    Extract name and contact from owner
    :param this_sample_id:
    :param owner:
    :param verbose: more output
    :return:
    """

    name = ""
    email = ""

    for child in owner:
        if "Name" == child.tag and child.text:
            name = child.text
        elif "Contact" == child.tag and "email" in child.attrib:
            email = child.attrib['email']

    return name, email


def parse_status(this_sample_id, status, verbose=False):
    """
    What is the status of this
    :param this_sample_id:
    :param status:
    :param verbose: more output
    :return:
    """

    if 'status' in status.attrib and status.attrib['status'] != 'live':
        sys.stderr.write("WARNING: {}. Status is not live: status: {}\n".format(this_sample_id, status.attrib['status']))
    if 'when' not in status.attrib:
        sys.stderr.write("WARNING: {}. No time stamp for status.\n".format(this_sample_id))
        return "unknown"
    return status.attrib['when']

def parse_biosample(biosample, verbose=False):
    """
    Parse a biosample object and print the information
    :param biosample: the biosample xml to parse
    :param verbose: more output
    :return:
    """

    this_sample_id = "Unknown"
    if 'accession' in biosample.attrib:
        this_sample_id = biosample.attrib['accession']
    else:
        sys.stderr.write(f"Unknown accession number in {biosample}, skipped\n")
        return None

    if verbose:
        message(f"Parsing {this_sample_id}", 'GREEN')

    global known_attrs_set
    global known_titles_set
    global known_attributes_set

    for attr in biosample.attrib:
        if (attr not in known_attrs_set):
            known_attrs_set.add(attr)

    contents = {x:"" for x in known_titles}

    # parse the known_attributes
    for x in known_attrs:
        if x in biosample.attrib:
            contents[x]=biosample.attrib[x]
        else:
            contents[x] = None
            sys.stderr.write(f"No attribute {x} in biosample\n")

    # parse the known_titles elements
    attributes = {}
    for child in biosample:
        if 'Ids' == child.tag:
            contents['Ids'] = parse_ids(this_sample_id, child)
        elif 'Description' == child.tag:
            (contents['Description - Title'], contents['Description - Comment'], contents['Organism']) = parse_description(this_sample_id, child)
        elif 'Attributes' == child.tag:
            for attr in child:
                (n,v) = parse_attribute(this_sample_id, attr)
                if not v:
                    sys.stderr.write(f"ERROR: no value returned for {n} from {this_sample_id}\n")
                    continue
                if n not in known_attributes_set:
                    known_attributes_set.add(n)
                if n in attributes:
                    if attributes[n] != v:
                        # sys.stderr.write("Redundant Atributes: ({}): Appending attribute value for {}. Previously had {} and now have {}\n".format(this_sample_id, n, attributes[n], v))
                        attributes[n] = attributes[n] + "; " + v.strip()
                else:
                    attributes[n]=v.strip()
        elif 'Owner' == child.tag:
            contents['Owner - Name'], contents['Owner - Email'] = parse_owner(this_sample_id, child)
        elif 'Status' == child.tag:
            contents['Release Date'] = parse_status(this_sample_id, child)
        elif 'Links' == child.tag:
            contents['Links'] = parse_links(this_sample_id, child)
        elif 'Models' == child.tag or 'Package' == child.tag:
            # this doesn't seem to contain a lot of information?
            continue
        else:
            sys.stderr.write("Skipped the title: {} as we don't know how to parse it\n".format(child.tag))

    
        contents['attributes'] = attributes

    return this_sample_id, contents



def write_outputs(data, minocc=1, verbose=False):
    """
    Write all the data out
    :param data: dict with keys = sample id, values = contents
    :param minocc: minimum occurrence of any key in the contents
    :param verbose: more output
    :return: nothing
    """

    # figure out the contents
    kcounts = {}
    acounts = {}
    for s in data:
        for k in data[s]:
            kcounts[k] = kcounts.get(k, 0)+1
        for k in data[s]['attributes']:
            acounts[k] = acounts.get(k, 0)+1

    kt = sorted(list(filter(lambda x: kcounts[x] > minocc, known_titles_set)))
    ka = sorted(list(filter(lambda x: acounts[x] > minocc, known_attributes_set)))


    sys.stdout.write("Sample Accession\t")
    sys.stdout.write("\t".join(kt))
    sys.stdout.write("\t")
    sys.stdout.write("\t".join(ka))
    sys.stdout.write("\n")


    for s in data:
        sys.stdout.write(s)
        for k in kt:
            if k in data[s]:
                sys.stdout.write(f"\t{data[s][k]}")
            else:
                sys.stdout.write("\t")
        for k in ka:
            if k in data[s]['attributes']:
                sys.stdout.write(f"\t{data[s]['attributes'][k]}")
            else:
                sys.stdout.write("\t")
        sys.stdout.write("\n")



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse an xml file of metadata for crAssphage')
    parser.add_argument('-d', help='directory of xml files to parse (e.g. /data/SRA/biosample/crassphage.xml)', required=True)
    parser.add_argument('-m', help='minimum number of files an attribute has to appear in [Default %(default)d]', type=int, default=1)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    data = {}
    for f in os.listdir(args.d):
        tree = ET.parse(os.path.join(args.d, f))
        biosampleset = tree.getroot()
        for biosample in biosampleset:
            tid, cont = parse_biosample(biosample)
            data[tid] = cont

    write_outputs(data, args.m)

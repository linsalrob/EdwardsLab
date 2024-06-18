"""
Parse an ENA XML file and print the data as a list
"""

import os
import sys
import argparse

__author__ = 'Rob Edwards'

import xml.etree.ElementTree as ET

# Function to extract data from XML and return as a list
def extract_data_from_xml(xml_content):
    root = ET.fromstring(xml_content)

    study = root.find('STUDY')
    accession = study.get('accession')
    alias = study.get('alias')
    center_name = study.get('center_name')
    broker_name = study.get('broker_name')

    identifiers = study.find('IDENTIFIERS')
    primary_id = identifiers.find('PRIMARY_ID').text
    secondary_id = identifiers.find('SECONDARY_ID').text
    submitter_id = identifiers.find('SUBMITTER_ID').text

    descriptor = study.find('DESCRIPTOR')
    study_title = descriptor.find('STUDY_TITLE').text
    study_abstract = descriptor.find('STUDY_ABSTRACT').text
    study_description = descriptor.find('STUDY_DESCRIPTION').text
    center_project_name = descriptor.find('CENTER_PROJECT_NAME').text
    study_type = descriptor.find('STUDY_TYPE').get('existing_study_type')

    study_links = study.find('STUDY_LINKS')
    study_link_list = []
    for link in study_links.findall('STUDY_LINK'):
        db = link.find('XREF_LINK/DB').text
        id = link.find('XREF_LINK/ID').text
        study_link_list.append(f"{db}:{id}")

    study_attributes = study.find('STUDY_ATTRIBUTES')
    study_attribute_list = []
    for attribute in study_attributes.findall('STUDY_ATTRIBUTE'):
        tag = attribute.find('TAG').text
        value = attribute.find('VALUE').text
        study_attribute_list.append(f"{tag}:{value}")

    # Combine all data into a single list
    data = [
        accession, alias, center_name, broker_name, primary_id,
        secondary_id, submitter_id, study_title, study_abstract,
        study_description, center_project_name, study_type,
        ";".join(study_link_list), ";".join(study_attribute_list)
    ]
    return data

def print_header():
    header = [
        "Accession", "Alias", "Center Name", "Broker Name", "Primary ID",
        "Secondary ID", "Submitter ID", "Study Title", "Study Abstract",
        "Study Description", "Center Project Name", "Study Type",
        "Study Links", "Study Attributes"
    ]

    print("\t".join(header))


def print_data(data):
    print("\t".join(data))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-d', help='directory of xml files', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    print_header()
    for xmlfile in os.listdir(args.d):
        if xmlfile.endswith(".xml"):
            if args.v:
                print("Parsing {xmlfile", file=sys.stderr)
            with open(os.path.join(args.d, xmlfile), 'r') as f:
                xml_content = f.read()
                data = extract_data_from_xml(xml_content)
                print_data(data)
        else:
            print("Skipping {xmlfile} because it does not end xml", file=sys.stderr)

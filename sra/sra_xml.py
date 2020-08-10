"""
Parse the SRA XML and extract useful information!


These are all the tags (in the current crassphage.xml file. Use sra_xml_print_all_attributes to get the list from an xml file):
Attribute	attribute_name
Attribute	display_name
Attribute	harmonized_name
Attributes
Attribute	unit
BioSample	access
BioSample	accession
BioSample	id
BioSample	last_update
BioSample	publication_date
BioSample	submission_date
Comment
Contact	email
Contact	lab
Contact	phone
Contacts
Description
First
Id	db
Id	db_label
Id	is_hidden
Id	is_primary
Ids
Last
Link	label
Links
Link	target
Link	type
Middle
Model
Models
Model	version
Name
Name	abbreviation
Name	url
OrganismName
Organism	taxonomy_id
Organism	taxonomy_name
Owner
Package	display_name
Paragraph
Status	status
Status	when
Title



"""

import os
import sys
import argparse
import xml.etree.ElementTree as ET

# we make a list so that the order is guaranteed, and then make a set for O(1) lookup
known_attrs = ['id', 'accession', 'last_update', 'access', 'publication_date', 'submission_date']
known_titles = ['Ids', 'Description - Title', 'Description - Comment', 'Owner - Name', 'Owner - Email', 'Release Date',
                'Links']
known_attributes = ['abs_air_humidity', 'age', 'age range', 'air_temp', 'alias', 'alkalinity', 'alt_elev', 'alternate_id', 'altitude', 'analyte_type', 'anonymized name', 'anonymized_name', 'assembly', 'assembly_method', 'assembly_method_and_version', 'assembly_method_version', 'assemblyversion', 'attribute_package', 'biochem_oxygen_dem', 'biomaterial_provider', 'biospecimen_repository', 'biospecimen_repository_sample_id', 'biotic_relationship', 'biotype', 'body sample site', 'body sample subsite', 'body_habitat', 'body_mass_index', 'body_product', 'breed', 'broker name', 'build_occup_type', 'building_setting', 'carb_dioxide', 'cell shape', 'chem_administration', 'chem_oxygen_dem', 'clone', 'collected_by', 'collection_date', 'common name', 'completeness score', 'completeness_estimated', 'contamination score', 'contamination_estimated', 'cultivar', 'culture_collecction', 'culture_collection', 'depth', 'derived_from', 'description', 'dev_stage', 'diet', 'disease', 'disease_stage', 'dreived-from', 'elev', 'ena checklist', 'ena first public', 'ena last update', 'ena-checklist', 'ena-first-public', 'ena-last-update', 'encoded_traits', 'env_biome', 'env_broad_scale', 'env_feature', 'env_local_scale', 'env_material', 'env_medium', 'env_package', 'environment', 'environmental-sample', 'environmental_sample', 'estimated_size', 'ethnicity', 'external id', 'extrachrom_elements', 'family_relationship', 'filter_type', 'finishing strategy (depth of coverage)', 'foodon ontology term', 'funding program', 'gap_accession', 'gap_consent_code', 'gap_consent_short_name', 'gap_sample_id', 'gap_subject_id', 'gastrointest_disord', 'gene calling method', 'genome size (bp)', 'genome_coverage', 'genotype', 'genus', 'geo_loc_name', 'geographic location (region and locality)', 'gold stamp id', 'gram staining', 'greengenes id', 'habitat', 'health_state', 'heat_cool_type', 'hhs_region', 'history of isolate', 'host', 'host_age', 'host_body_habitat', 'host_body_mass_index', 'host_body_product', 'host_body_temp', 'host_description', 'host_diet', 'host_disease', 'host_disease_outcome', 'host_disease_stage', 'host_family_relationship', 'host_genotype', 'host_health_state', 'host_height', 'host_last_meal', 'host_life_stage', 'host_occupation', 'host_phenotype', 'host_pulse', 'host_sex', 'host_subject_id', 'host_taxid', 'host_tissue_sampled', 'host_tot_mass', 'identification method', 'identified_by', 'ifsac+ category', 'ihmc_medication_code', 'indoor_space', 'insdc center alias', 'insdc center name', 'insdc first public', 'insdc last update', 'insdc status', 'investigation_type', 'isol_growth_condt', 'isolate', 'isolate_name_alias', 'isolation comments', 'isolation site', 'isolation_source', 'lab_host', 'label', 'lat_lon', 'light_type', 'liver_disord', 'locus_tag_prefix', 'mapping_method', 'mapping_method_and_version', 'mapping_method_version', 'mating_type', 'medic_hist_perform', 'metagenome-source', 'metagenome_source', 'metagenomic', 'misc_param', 'misc_param: hmp body site', 'misc_param: hmp supersite', 'molecular_data_type', 'motility', 'nitrate', 'note', 'nucleic acid extraction', 'num_replicons', 'number_of_identified_antimicrobial_resistance_genes', 'occup_samp', 'occupant_dens_samp', 'organism_count', 'orgmod_note', 'oxy_stat_samp', 'oxygen requirement', 'passage_history', 'pathogenicity', 'pathotype', 'pathovar', 'perturbation', 'pfge_primaryenzyme_pattern', 'pfge_secondaryenzyme_pattern', 'ph', 'phenotypes', 'phosphate', 'pre_treatment', 'project_name', 'project_type', 'projectaccession', 'propagation', 'publicaccession', 'quality_assessment_method', 'quality_assessment_method_and_version', 'quality_assessment_method_version', 'race', 'reactor_type', 'ref_biomaterial', 'region (hhs)', 'rel_air_humidity', 'rel_to_oxygen', 'salinity', 'samp_collect_device', 'samp_mat_process', 'samp_salinity', 'samp_size', 'samp_store_dur', 'samp_store_loc', 'samp_store_temp', 'samp_vol_we_dna_ext', 'sample comment', 'sample derived from', 'sample number', 'sample_identifier', 'sample_name', 'sample_type', 'sequencing depth', 'sequencing method', 'sequencing_meth', 'sequencingtechnology', 'serogroup', 'serotype', 'serovar', 'sewage_type', 'sex', 'sludge_retent_time', 'smoker', 'sop', 'source_material_id', 'space_typ_state', 'special_diet', 'species', 'specimen_voucher', 'sporulation', 'sra accession', 'store_cond', 'strain', 'strain_name_alias', 'study_design', 'study_disease', 'study_name', 'sub_species', 'subgroup', 'subject_id', 'subject_is_affected', 'submitted_sample_id', 'submitted_subject_id', 'submitter id', 'submitter_handle', 'subsrc_note', 'subtype', 'supplier_name', 'suspend_solids', 'temp', 'temperature optimum', 'temperature range', 'timepoint', 'tissue', 'title', 'tot_phosphate', 'treatment', 'trophic_level', 'typ_occupant_dens', 'type', 'type-material', 'type_strain', 'ukzn', 'uploaddate', 'value', 'ventilation_type', 'wastewater_type']

known_attrs_set = set(known_attrs)
known_titles_set = set(known_titles)
known_attributes_set = set(known_attributes)


def parse_ids(this_sample_id, ids):
    """
    Parse the IDs field
    :param this_sample_id: The current sample ID for debugging
    :param ids: The ID xml
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
    return id_text


def parse_links(this_sample_id, links):
    """
    Parse the links field
    :param this_sample_id:
    :param links:
    :return:
    """

    link_text=""
    for child in links:
        if 'type' in child.attrib and "entrez" == child.attrib['type']:
            link_text = link_text + "{}|{}; ".format(child.attrib['target'], child.text)
        elif 'type' in child.attrib and "url" == child.attrib['type']:
            link_text = link_text + "url|[{}]({}); ".format(child.text, child.attrib['label'])

    return link_text

def parse_description(this_sample_id, desc):
    """
    Parse the description field and return the title and the
    :param this_sample_id: The current sample ID for debugging
    :param desc: The description xml
    :return:
    """
    desc_title = "";
    desc_comment = "";
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

    return (desc_title, desc_comment, desc_org)


def parse_attribute(this_sample_id, attr):
    """
    Parse out the attribute field
    :param this_sample_id:
    :param attr:
    :return:
    """

    if 'harmonized_name' in attr.attrib:
        return attr.attrib['harmonized_name'], attr.text
    elif 'attribute_name' in attr.attrib:
        return attr.attrib['attribute_name'].lower(), attr.text
    else:
        sys.stderr.write("WARNING: ({}): Neither harmonized_name nor attribute name in attribute {}".format(this_sample_id, attr.text))
        return "",""


def parse_owner(this_sample_id, owner):
    """
    Extract name and contact from owner
    :param this_sample_id:
    :param owner:
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


def parse_status(this_sample_id, status):
    """
    What is the status of this
    :param this_sample_id:
    :param status:
    :return:
    """

    if 'status' in status.attrib and status.attrib['status'] != 'live':
        sys.stderr.write("WARNING: {}. Status is not live: status: {}\n".format(this_sample_id, status.attrib['status']))
    if 'when' not in status.attrib:
        sys.stderr.write("WARNING: {}. No time stamp for status.\n".format(this_sample_id))
        return "unknown"
    return status.attrib['when']

def parse_biosample(biosample, header):
    """
    Parse a biosample object and print the information
    :param biosample: the biosample xml to parse
    :param header: an integer. We only print the header line when header == 0
    :return:
    """
    # this is used for debugging
    this_sample_id = "Unknown"
    if 'accession' in biosample.attrib:
        this_sample_id = biosample.attrib['accession']
    else:
        sys.stderr.write(f"Unknown accession number in {biosample}\n")

    global known_attrs_set
    global known_titles_set
    global known_attributes_set


    for attr in biosample.attrib:
        if (attr not in known_attrs_set):
            sys.stderr.write("WARNING: ({}) New key in biosample: {}. Please add this to the known_attrs variable\n".format(this_sample_id, attr))
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
                if n not in known_attributes_set:
                    # sys.stderr.write("ERROR: New attribute {} in {}. Please append to `known_attributes` in this code\n".format(n, this_sample_id))
                    known_attributes_set.add(n)
                if n in attributes:
                    if attributes[n] != v:
                        sys.stderr.write("Redundant Atributes: ({}): Appending attribute value for {}. Previously had {} and now have {}\n".format(this_sample_id, n, attributes[n], v))
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

    if header == 0:
        sys.stdout.write("Sample Accession\t")
        sys.stdout.write("\t".join(known_titles))
        sys.stdout.write("\t")
        sys.stdout.write("\t".join(known_attributes))
        sys.stdout.write("\n")


    sys.stdout.write(this_sample_id)
    for k in known_titles:
        sys.stdout.write("\t{}".format(contents[k]))
    for a in known_attributes:
        if a in attributes:
            sys.stdout.write("\t{}".format(attributes[a]))
        else:
            sys.stdout.write("\t")
    sys.stdout.write("\n")



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse an xml file of metadata for crAssphage')
    parser.add_argument('-f', help='xml file to parse (e.g. /data/SRA/biosample/crassphage.xml)', required=True)
    args = parser.parse_args()

    tree = ET.parse(args.f)
    biosampleset = tree.getroot()
    header=0
    for biosample in biosampleset:
        parse_biosample(biosample, header)
        header=header+1

    if len(known_attrs_set) != len(known_attrs):
        sys.stderr.write("PLEASE UPDATE THIS CODE. Edit known_attrs with this list:\n")
        l = sorted(list(known_attrs_set))
        sys.stderr.write(f"{l}\n")
    if len(known_attributes_set) != len(known_attributes):
        sys.stderr.write("PLEASE UPDATE THIS CODE. Edit known_attributes with this list:\n")
        l = sorted(list(known_attributes_set))
        sys.stderr.write(f"s.update({l})\n")


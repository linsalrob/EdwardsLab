"""
Create the NCBI BioSample submission file from our spreadsheet file
"""

import os
import sys
import argparse
import re
from roblib import bcolors

countries = {"Afghanistan", "Albania", "Algeria", "American Samoa", "Andorra", "Angola", "Anguilla", "Antarctica", "Antigua and Barbuda", "Arctic Ocean", "Argentina", "Armenia", "Aruba", "Ashmore and Cartier Islands", "Atlantic Ocean", "Australia", "Austria", "Azerbaijan", "Bahamas", "Bahrain", "Baltic Sea", "Baker Island", "Bangladesh", "Barbados", "Bassas da India", "Belarus", "Belgium", "Belize", "Benin", "Bermuda", "Bhutan", "Bolivia", "Borneo", "Bosnia and Herzegovina", "Botswana", "Bouvet Island", "Brazil", "British Virgin Islands", "Brunei", "Bulgaria", "Burkina Faso", "Burundi", "Cambodia", "Cameroon", "Canada", "Cape Verde", "Cayman Islands", "Central African Republic", "Chad", "Chile", "China", "Christmas Island", "Clipperton Island", "Cocos Islands", "Colombia", "Comoros", "Cook Islands", "Coral Sea Islands", "Costa Rica", "Cote d'Ivoire", "Croatia", "Cuba", "Curacao", "Cyprus", "Czech Republic", "Democratic Republic of the Congo", "Denmark", "Djibouti", "Dominica", "Dominican Republic", "East Timor", "Ecuador", "Egypt", "El Salvador", "Equatorial Guinea", "Eritrea", "Estonia", "Ethiopia", "Europa Island", "Falkland Islands (Islas Malvinas)", "Faroe Islands", "Fiji", "Finland", "France", "French Guiana", "French Polynesia", "French Southern and Antarctic Lands", "Gabon", "Gambia", "Gaza Strip", "Georgia", "Germany", "Ghana", "Gibraltar", "Glorioso Islands", "Greece", "Greenland", "Grenada", "Guadeloupe", "Guam", "Guatemala", "Guernsey", "Guinea", "Guinea-Bissau", "Guyana", "Haiti", "Heard Island and McDonald Islands", "Honduras", "Hong Kong", "Howland Island", "Hungary", "Iceland", "India", "Indian Ocean", "Indonesia", "Iran", "Iraq", "Ireland", "Isle of Man", "Israel", "Italy", "Jamaica", "Jan Mayen", "Japan", "Jarvis Island", "Jersey", "Johnston Atoll", "Jordan", "Juan de Nova Island", "Kazakhstan", "Kenya", "Kerguelen Archipelago", "Kingman Reef", "Kiribati", "Kosovo", "Kuwait", "Kyrgyzstan", "Laos", "Latvia", "Lebanon", "Lesotho", "Liberia", "Libya", "Liechtenstein", "Line Islands", "Lithuania", "Luxembourg", "Macau", "Madagascar", "Malawi", "Malaysia", "Maldives", "Mali", "Malta", "Marshall Islands", "Martinique", "Mauritania", "Mauritius", "Mayotte", "Mediterranean Sea", "Mexico", "Micronesia", "Midway Islands", "Moldova", "Monaco", "Mongolia", "Montenegro", "Montserrat", "Morocco", "Mozambique", "Myanmar", "Namibia", "Nauru", "Navassa Island", "Nepal", "Netherlands", "New Caledonia", "New Zealand", "Nicaragua", "Niger", "Nigeria", "Niue", "Norfolk Island", "North Korea", "North Sea", "Northern Mariana Islands", "Norway", "Oman", "Pacific Ocean", "Pakistan", "Palau", "Palmyra Atoll", "Panama", "Papua New Guinea", "Paracel Islands", "Paraguay", "Peru", "Philippines", "Pitcairn Islands", "Poland", "Portugal", "Puerto Rico", "Qatar", "Republic of the Congo", "Reunion", "Romania", "Ross Sea", "Russia", "Rwanda", "Saint Helena", "Saint Kitts and Nevis", "Saint Lucia", "Saint Pierre and Miquelon", "Saint Vincent and the Grenadines", "Samoa", "San Marino", "Sao Tome and Principe", "Saudi Arabia", "Senegal", "Serbia", "Seychelles", "Sierra Leone", "Singapore", "Sint Maarten", "Slovakia", "Slovenia", "Solomon Islands", "Somalia", "South Africa", "South Georgia and the South Sandwich Islands", "South Korea", "South Sudan", "Southern Ocean", "Spain", "Spratly Islands", "Sri Lanka", "State of Palestine", "Sudan", "Suriname", "Svalbard", "Swaziland", "Sweden", "Switzerland", "Syria", "Taiwan", "Tajikistan", "Tanzania", "Tasman Sea", "Thailand", "The former Yugoslav Republic of Macedonia", "Togo", "Tokelau", "Tonga", "Trinidad and Tobago", "Tromelin Island", "Tunisia", "Turkey", "Turkmenistan", "Turks and Caicos Islands", "Tuvalu", "USA", "Uganda", "Ukraine", "United Arab Emirates", "United Kingdom", "Uruguay", "Uzbekistan", "Vanuatu", "Venezuela", "Viet Nam", "Virgin Islands", "Wake Island", "Wallis and Futuna", "West Bank", "Western Sahara", "Yemen", "Zambia", "Zimbabwe"}

# "
# Add terms that identify the major environment type(s) where your sample was collected. Recommend subclasses of biome
# [ENVO:00000428]. Multiple terms can be separated by one or more pipes
# e.g.:  mangrove biome [ENVO:01000181]|estuarine biome [ENVO:01000020] 
env_broad_scale = {
    "WWTP": "ENVO:00002001",
    'fecal sample': 'ENVO:01001029',
    'post-STP, refugee camp': 'ENVO:00002001|NCIT:C85867',
    'pre-sewage treatment plant' : 'ENVO:00002001',
    'Raw Sewage' : 'ENVO:00002001'
}

# Add terms that identify environmental entities having causal influences upon the entity at time of sampling,
# multiple terms can be separated by pipes, e.g.:  shoreline [ENVO:00000486]|intertidal zone [ENVO:00000316]
env_local_scale = {
    "WWTP": "ENVO:00002018",
    'fecal sample' : 'ENVO:00002003',
    'post-STP, refugee camp': 'ENVO:00002018',
    'pre-sewage treatment plant' : 'ENVO:00002018',
    'Raw Sewage' : 'ENVO:00002018'
}

# Add terms that identify the material displaced by the entity at time of sampling. Recommend subclasses of
# environmental material [ENVO:00010483]. Multiple terms can be separated by
# pipes e.g.: estuarine water [ENVO:01000301]|estuarine mud [ENVO:00002160]
env_medium = {
    "WWTP": "ENVO:00003043",
    'fecal sample' : 'UBERON:0001988',
    'post-STP, refugee camp': 'ENVO:00003043',
    'pre-sewage treatment plant' : 'ENVO:00003043',
    'Raw Sewage' : 'ENVO:00003043'
}

predefined = {
    "SRA" : {"bioproject_accession" : "PRJNA510571", "organism" : "uncultured crAssphage"},
    "WWTP" : {"bioproject_accession" : "PRJNA510571", "organism" : "uncultured crAssphage", "wastewater_type" : "human waste", "sewage_type" : "municiple"},
    'fecal sample': {"bioproject_accession" : "PRJNA510571", "organism" : "uncultured crAssphage"},
    'post-STP, refugee camp': {"bioproject_accession" : "PRJNA510571", "organism" : "uncultured crAssphage", "wastewater_type" : "human waste", "sewage_type" : "municiple"},
    'pre-sewage treatment plant': {"bioproject_accession" : "PRJNA510571", "organism" : "uncultured crAssphage", "wastewater_type" : "human waste", "sewage_type" : "municiple"},
    'Raw Sewage': {"bioproject_accession" : "PRJNA510571", "organism" : "uncultured crAssphage", "wastewater_type" : "human waste", "sewage_type" : "municiple"}
}

columns = {
    "date" : "collection_date", "latlon" : "lat_lon", "name" : "sample_title",
    "address" : "address", "altitude" : "altitude", "method" : "extraction_method", "sex" : "host_sex",
    "volunteer" : "host_subject_id", "note" : "note", "contact" : "provider",
    "sampletype" : "sample_frequency", "description" : "description",
    "university" : "university"
}

def check_existing(data, tid, tag, newval):
    """
    Check the new value compared to an existing value
    :param data: the data dictionary
    :param tid: the id of this object
    :param tag: the column header
    :param newval: the new value
    :return:
    """

    if tid not in data:
        sys.stderr.write(f"{bcolors.FAIL}FATAL: Trying to add {tid} and {tag} but don't have {tid} yet!{bcolors.ENDC}\n")
        sys.exit(-1)

    if tag not in data[tid]:
        return

    if data[tid][tag] == newval:
        return

    sys.stderr.write(f"{bcolors.WARNING}WARNING: for {tid} tag: {tag} we previously had {data[tid][tag]} but are replacing it with {newval}!{bcolors.ENDC}\n")


def read_attributes(filename = "attributes.txt"):
    """
    Read the SRA Attributes
    """
        
    attributes = {}
    with open(filename, 'r', encoding='utf-8') as f:
        for l in f:
            if l.endswith("\t\n"):
                # no values for this attributes
                attributes[l.strip()] = {}
                continue
            (sid, attr) = l.strip().split("\t")
            attributes[sid] = {}
            for p in attr.split("||"):
                #sys.stderr.write(p)
                if ":" not in p:
                    continue
                if 2 == p.count(":"):
                    (k, v, units) = p.split(":")
                elif 1 == p.count(":"):
                    (k, v) = p.split(":")
                else:
                    sys.stderr.write("Error unpacking |{}|\n".format(p))
                    continue
                k = k.strip()
                v = v.strip()
                k = k.replace(" ", "_")
                k = k.replace("-", "_")
                k = k.lower()
                # correct some wrong attributes
                if k in ['country', 'geo loc name', 'geographic location (country and/or sea region)', 
                         'geographic location (country and/or sea)', 'geographic location (country and/or sea, region)',
                         'geographic location (country and/or sea,region)', 'geographic location (country)', 
                         'geographic location (country:region,area)', 'geographic location (locality)', 
                         'geographic locations', 'geographic origin', 'geographical location', 
                         'geographical location (country:region, location)', 'geolocname']:
                    k = 'geo_loc_name'
                if k == 'sex' or k == 'gender':
                    k = 'host_sex'
                if k in {'/environmental_sample', 'geographic_location_(country'}:
                    continue
                if 'geographic_location' in k:
                    continue
                attributes[sid][k] = v


    return attributes

def parse_file(filename, attributes):
    """
    PArse a tsv file
    :param filename: what is the name
    :param attributes: what are the sra attributes
    :return:
    """

    data = {}

    with open(filename, 'r') as f:
        headers = []
        sourceidx = -1
        for l in f:
            p = l.strip().split("\t")
            if l.startswith("sequence ID"):
                if "source" not in p:
                    sys.stderr.write(f"{bcolors.FAIL}FATAL:{bcolors.ENDC} No source found in {f}\n")
                    sys.exit(-1)
                sourceidx = p.index("source")
                headers = p
                continue

            if p[0] not in data:
                data[p[0]] = {
                    "bioproject_accession" : "PRJNA510571",
                    "organism" : "uncultured crAssphage",
                    "cultivar" : "not applicable",
                    "ecotype" : "not applicable", 
                    "isolate" : "not applicable",
                    "strain" : "not applicable",
                    "isolation and growth condition" : "not applicable"

                }

            # deal with the country
            if "country" not in headers:
                sys.stderr.write(f"{bcolors.FAIL}FATAL:{bcolors.ENDC} No country found in {f}\n")
                sys.exit(-1)
            cidx = headers.index("country")
            if p[cidx] and p[cidx] not in countries:
                sys.stderr.write(f"{bcolors.FAIL}FATAL:{bcolors.ENDC}: |{p[cidx]}| is not a valid country\n")
                sys.exit(-1)
            if p[cidx]:
                mycountry = p[cidx]

                # now add the locality to this
                if "locality" not in headers:
                    sys.stderr.write(f"{bcolors.FAIL}FATAL:{bcolors.ENDC} No locality found in {f}\n")
                    sys.exit(-1)
                lidx = headers.index("locality")
                if p[lidx]:
                    mycountry += f": {p[lidx]}"
                    
                check_existing(data, p[0], 'geo_loc_name', mycountry)
                data[p[0]]['geo_loc_name'] = mycountry


            # deal with the source
            src = p[sourceidx]
            if "SRA" == src:
                sra_id = p[headers.index('name')].split(" ")[2]
                if not re.match('^[SDE]RR\d+$', sra_id):
                    sys.stderr.write(f"{bcolors.FAIL}SRA ID{bcolors.ENDC}: {sra_id}\n")
                # Uncultivated Euryarchaeota archaeon UBA41 genome recovered from ERX556009
                sn = f"Uncultured crAssphage amplicon recovered from {sra_id}"
                check_existing(data, p[0], 'sample_name', sn)
                data[p[0]]['sample_name'] = sn
                p[sourceidx] = sn
                data[p[0]]['extraction_method'] = f"Sequence was extracted from bioproject {sra_id} using Gretel version 0.08."
                data[p[0]]['sra_id'] = sra_id
            else:
                sn = f"Uncultured crAssphage amplicon"
                contaxtidx = headers.index('contact')
                if p[contaxtidx]:
                    sn += f" provided by {p[contaxtidx]}"
                    #uidx = headers.index('university')
                    #if p[uidx]:
                    #    sn += f" at {p[uidx]}"
                nameidx = headers.index('name')
                p[nameidx] = p[nameidx].replace('.', ' ')
                sn += f". Sample {p[nameidx]}"
                if len(sn) > 128:
                    sys.stderr.write("f{bcolor.RED}FATAL too long{bcolor.WHITE} {sn}\n")
                data[p[0]]['sample'] = p[nameidx]
                sn = re.sub('primer\s*[abc]', '', sn, flags=re.I)
                sn = sn.replace('  ', ' ')
                sn = sn.replace('_', ' ')
                check_existing(data, p[0], 'sample_name', sn)
                data[p[0]]['sample_name'] = sn
                # deal with envO terms
                if src not in env_broad_scale:
                    sys.stderr.write(f"{bcolors.FAIL}FATAL:{bcolors.ENDC}: |{src}| is not a valid source\n")
                    sys.exit(-1)
                check_existing(data, p[0], 'env_broad_scale', env_broad_scale[src])
                data[p[0]]['env_broad_scale'] = env_broad_scale[src]
                check_existing(data, p[0], 'env_local_scale', env_local_scale[src])
                data[p[0]]['env_local_scale'] = env_local_scale[src]
                check_existing(data, p[0], 'env_medium', env_medium[src])
                data[p[0]]['env_medium'] = env_medium[src]

                # predefined terms
                for pd in predefined[src]:
                    check_existing(data, p[0], pd, predefined[src][pd])
                    data[p[0]][pd] = predefined[src][pd]

            # now the rest:
            for c in columns:
                if c not in headers:
                    sys.stderr.write(f"{bcolors.WARNING}WARN:{bcolors.ENDC} {c} not found in {f}\n")
                    continue
                if "method" == c and "SRA" == src:
                    # we set this above and add the SRA ID to it
                    continue
                idx = headers.index(c)
                if p[idx]:
                    check_existing(data, p[0], columns[c], p[idx])
                    data[p[0]][columns[c]] = p[idx]
            data[p[0]]['src'] = src

            # finally, if it is an SRA id, we trust their attributes better than ours :)
            if 'sra_id' in data[p[0]] and data[p[0]]['sra_id'] in attributes:
                data[p[0]].update(attributes[data[p[0]]['sra_id']])

    return data

def fix_data(data, verbose=False):
    """
    Fix some of the issues with the data
    """
    
    """
    columns = {
        "date" : "collection_date", "latlon" : "lat_lon", "name" : "sample_title",
        "address" : "address", "altitude" : "altitude", "method" : "extraction_method", "sex" : "host_sex",
        "volunteer" : "host_subject_id", "note" : "note", "contact" : "provider",
        "sampletype" : "sample_frequency", "description" : "description",
        "university" : "university"
    }
    """

    for d in data:
        #if "Ghan_20120101_Australia_10" == d:
        #    sys.stderr.write("Yes, Ghan_20120101_Australia_10 really here\n\n\n")
        #    sys.stderr.write("KEYS: {}\n".format("\n".join(data[d].keys())))
        #    sys.stderr.write("Yes, Ghan_20120101_Australia_10 really here\n\n\n")
        if "host_sex" in data[d] and "m" == data[d]['host_sex'].lower():
            data[d]['host_sex'] = "male"

        if "host_sex" in data[d] and "f" == data[d]['host_sex'].lower():
            data[d]['host_sex'] = "female"
        
        if 'collection_date' in data[d]:
            if 'missing' == data[d]['collection_date']:
                pass
            elif 10 == len(data[d]['collection_date']) and 2 == data[d]['collection_date'].count('-'):
                pass
            elif 9 == len(data[d]['collection_date']) and (1 == data[d]['collection_date'].count('/') or 1 == data[d]['collection_date'].count('-')):
                pass
            elif 4 == len(data[d]['collection_date']):
                pass
            elif 8 == len(data[d]['collection_date']):
                data[d]['collection_date'] = data[d]['collection_date'][0:4] + "-" + data[d]['collection_date'][4:6] + "-" + data[d]['collection_date'][6:8]
            else:
                #sys.stderr.write(f"{bcolors.WARNING}Date: {data[d]['collection_date']}{bcolors.ENDC}\n")
                pass
    
        for t in ['env_broad_scale', 'env_local_scale', 'env_medium', 'collection_date', 'lat_lon', 'geo_loc_name']:
            if t not in data[d] or not data[d][t]:
                data[d][t] = "missing"
        
    return data

def print_all_old(data, outputf, verbose=False):
    """
    Print out everything. This is the old version that does not remove empty columns

    Note that biosamples are limited to 1,000 entries per time, so we need to make
    three files.

    :param data: the dict of data
    :param outputf: the file to write
    :param verbose: more output
    :return:
    """

    # note we need to have a unique field that does not include
    # (excluding sample name, title, bioproject accession and description).

    order = ["sample_name", "bioproject_accession", "organism", "collection_date", \
            "env_broad_scale", "env_local_scale", "env_medium", "geo_loc_name", "sample", "host", "lat_lon"]

    # figure out all the keys
    ak = set()
    for d in data:
        ak.update(data[d].keys())
    [ak.remove(x) for x in order if x in ak]
    # not including a sample title.
    ak.remove('sample_title')
    # these are being stored for convenience so we can print them out
    ak.remove('src')
    ak.remove('sra_id')
    otherkeys = sorted(ak)


    seen = set()

    if '.tsv' in outputf:
        outputf = outputf.replace('.tsv', '')

    filecount = 0
    out = open(f"{outputf}.{filecount}.tsv", "w")
    out.write("{}\t{}\n".format("\t".join(order), "\t".join(otherkeys)))
    linecount = 0
    for d in data:
        if linecount >= 900:
            linecount = 0
            filecount += 1
            out.close()
            out = open(f"{outputf}.{filecount}.tsv", "w")
            out.write("{}\t{}\n".format("\t".join(order), "\t".join(otherkeys)))

        # out.write(d)
        if data[d]['sample_name'] in seen:
            continue
        data[d]['sample_title'] = re.sub('primer\s*[abc]', '', data[d]['sample_title'], flags=re.I)
        data[d]['sample_title'] = data[d]['sample_title'].replace('  ', ' ')
        for i, o in enumerate(order):
            if o not in data[d]:
                data[d][o] = ""
            if i > 0:
                out.write("\t")
            out.write(f"{data[d][o]}")
        for o in otherkeys:
            if o not in data[d]:
                data[d][o] = ""
            out.write(f"\t{data[d][o]}")
        out.write("\n")
        linecount += 1
        seen.add(data[d]['sample_name'])

    return data


def print_outputs(filecount, outputf, headers, outputs):
    """
    Write the output files
    :param filecount: the filenumber to write. We do not increment it here
    :param outputf: the output file base
    :param headers: the array of headers
    :param outputs: array of arrays
    :return: None
    """

    ## remove any empty columns
    counts = [0 for x in headers]
    for o in outputs:
        for i,j in enumerate(o):
            if j:
                counts[i] += 1

    with open(f"{outputf}.{filecount}.tsv", "w") as out:
        l = [headers[i] for i,j in enumerate(counts) if j > 0]
        out.write("{}\n".format("\t".join(l)))
        for o in outputs:
            l = [o[i] for i,j in enumerate(counts) if j > 0]
            out.write("{}\n".format("\t".join(l)))
    return


def print_all(data, outputf, verbose=False):
    """
    Print out everything. This version removes all the empty columns

    Note that biosamples are limited to 1,000 entries per time, so we need to make
    three files.

    :param data: the dict of data
    :param outputf: the file to write
    :param verbose: more output
    :return:
    """

    # note we need to have a unique field that does not include
    # (excluding sample name, title, bioproject accession and description).

    order = ["sample_name", "bioproject_accession", "organism", "collection_date", \
            "env_broad_scale", "env_local_scale", "env_medium", "geo_loc_name", "sample", "host", "lat_lon"]

    # figure out all the keys
    ak = set()
    for d in data:
        ak.update(data[d].keys())
    [ak.remove(x) for x in order if x in ak]
    # not including a sample title.
    ak.remove('sample_title')
    # these are being stored for convenience so we can print them out
    ak.remove('src')
    ak.remove('sra_id')
    otherkeys = sorted(ak)


    seen = set()

    if '.tsv' in outputf:
        outputf = outputf.replace('.tsv', '')

    filecount = 0

    headers = order + otherkeys
    outputs = []
    for d in data:
        if data[d]['sample_name'] in seen:
            continue
        data[d]['sample_title'] = re.sub('primer\s*[abc]', '', data[d]['sample_title'], flags=re.I)
        data[d]['sample_title'] = data[d]['sample_title'].replace('  ', ' ')
        thisoutput = []
        for i, o in enumerate(order):
            if o not in data[d]:
                data[d][o] = ""
            thisoutput.append(data[d][o])

        for o in otherkeys:
            if o not in data[d]:
                data[d][o] = ""
            thisoutput.append(data[d][o])
        seen.add(data[d]['sample_name'])
        outputs.append(thisoutput)
        if len(outputs) > 900:
            print_outputs(filecount, outputf, headers, outputs)
            filecount += 1
            outputs = []

    print_outputs(filecount, outputf, headers, outputs)

    return data

def print_sample_id_name(data, outputf, verbose=False):
    """
    Print out just the sample id and the name. this is so we can go back again later
    to connect biosample to sequences
    :param data: the dict of data
    :param outputf: the file to write
    :param verbose: more output
    :return:
    """
    with open(outputf, 'w') as out:
        for d in data:
            out.write("{}\t{}\n".format(d, data[d]['sample_name']))


__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse files for the NCBI BioSample Submission')
    parser.add_argument('-d', help='directory of csv files', required=True)
    parser.add_argument('-n', help='biosample file to write to submit to the NCBI', required=True)
    parser.add_argument('-a', help="sra attributes file, default=attributes.txt")
    parser.add_argument('-i', help='id/biosample name mapping file to go back to the sequences', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()
    
    attributes = {}
    if args.a:
        attributes = read_attributes(args.a)
    else:
        attributes = read_attributes()


    data = {}
    for f in os.listdir(args.d):
        sys.stderr.write(f"{bcolors.OKGREEN}PARSING:{bcolors.ENDC} {f}\n")
        newdata = parse_file(os.path.join(args.d, f), attributes)
        data.update(newdata)
        if 'DRR002659' in newdata:
            sys.stderr.write(f"{bcolors.WARNING}DRR002659 now in NEW data at main. Where did it come from?\n")
            sys.exit(-1)
        if 'DRR002659' in data:
            sys.stderr.write(f"{bcolors.WARNING}DRR002659 now in data at main. Where did it come from?\n")
            sys.exit(-1)


    data = fix_data(data, args.v)

    data = print_all(data, args.n, args.v)
    print_sample_id_name(data, args.i, args.v)

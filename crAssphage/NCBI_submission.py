"""
Create the NCBI BioSample submission file from our spreadsheet file
"""

import os
import sys
import argparse


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

countries = {"Afghanistan", "Albania", "Algeria", "American Samoa", "Andorra", "Angola", "Anguilla", "Antarctica", "Antigua and Barbuda", "Arctic Ocean", "Argentina", "Armenia", "Aruba", "Ashmore and Cartier Islands", "Atlantic Ocean", "Australia", "Austria", "Azerbaijan", "Bahamas", "Bahrain", "Baltic Sea", "Baker Island", "Bangladesh", "Barbados", "Bassas da India", "Belarus", "Belgium", "Belize", "Benin", "Bermuda", "Bhutan", "Bolivia", "Borneo", "Bosnia and Herzegovina", "Botswana", "Bouvet Island", "Brazil", "British Virgin Islands", "Brunei", "Bulgaria", "Burkina Faso", "Burundi", "Cambodia", "Cameroon", "Canada", "Cape Verde", "Cayman Islands", "Central African Republic", "Chad", "Chile", "China", "Christmas Island", "Clipperton Island", "Cocos Islands", "Colombia", "Comoros", "Cook Islands", "Coral Sea Islands", "Costa Rica", "Cote d'Ivoire", "Croatia", "Cuba", "Curacao", "Cyprus", "Czech Republic", "Democratic Republic of the Congo", "Denmark", "Djibouti", "Dominica", "Dominican Republic", "East Timor", "Ecuador", "Egypt", "El Salvador", "Equatorial Guinea", "Eritrea", "Estonia", "Ethiopia", "Europa Island", "Falkland Islands (Islas Malvinas)", "Faroe Islands", "Fiji", "Finland", "France", "French Guiana", "French Polynesia", "French Southern and Antarctic Lands", "Gabon", "Gambia", "Gaza Strip", "Georgia", "Germany", "Ghana", "Gibraltar", "Glorioso Islands", "Greece", "Greenland", "Grenada", "Guadeloupe", "Guam", "Guatemala", "Guernsey", "Guinea", "Guinea-Bissau", "Guyana", "Haiti", "Heard Island and McDonald Islands", "Honduras", "Hong Kong", "Howland Island", "Hungary", "Iceland", "India", "Indian Ocean", "Indonesia", "Iran", "Iraq", "Ireland", "Isle of Man", "Israel", "Italy", "Jamaica", "Jan Mayen", "Japan", "Jarvis Island", "Jersey", "Johnston Atoll", "Jordan", "Juan de Nova Island", "Kazakhstan", "Kenya", "Kerguelen Archipelago", "Kingman Reef", "Kiribati", "Kosovo", "Kuwait", "Kyrgyzstan", "Laos", "Latvia", "Lebanon", "Lesotho", "Liberia", "Libya", "Liechtenstein", "Line Islands", "Lithuania", "Luxembourg", "Macau", "Madagascar", "Malawi", "Malaysia", "Maldives", "Mali", "Malta", "Marshall Islands", "Martinique", "Mauritania", "Mauritius", "Mayotte", "Mediterranean Sea", "Mexico", "Micronesia", "Midway Islands", "Moldova", "Monaco", "Mongolia", "Montenegro", "Montserrat", "Morocco", "Mozambique", "Myanmar", "Namibia", "Nauru", "Navassa Island", "Nepal", "Netherlands", "New Caledonia", "New Zealand", "Nicaragua", "Niger", "Nigeria", "Niue", "Norfolk Island", "North Korea", "North Sea", "Northern Mariana Islands", "Norway", "Oman", "Pacific Ocean", "Pakistan", "Palau", "Palmyra Atoll", "Panama", "Papua New Guinea", "Paracel Islands", "Paraguay", "Peru", "Philippines", "Pitcairn Islands", "Poland", "Portugal", "Puerto Rico", "Qatar", "Republic of the Congo", "Reunion", "Romania", "Ross Sea", "Russia", "Rwanda", "Saint Helena", "Saint Kitts and Nevis", "Saint Lucia", "Saint Pierre and Miquelon", "Saint Vincent and the Grenadines", "Samoa", "San Marino", "Sao Tome and Principe", "Saudi Arabia", "Senegal", "Serbia", "Seychelles", "Sierra Leone", "Singapore", "Sint Maarten", "Slovakia", "Slovenia", "Solomon Islands", "Somalia", "South Africa", "South Georgia and the South Sandwich Islands", "South Korea", "South Sudan", "Southern Ocean", "Spain", "Spratly Islands", "Sri Lanka", "State of Palestine", "Sudan", "Suriname", "Svalbard", "Swaziland", "Sweden", "Switzerland", "Syria", "Taiwan", "Tajikistan", "Tanzania", "Tasman Sea", "Thailand", "The former Yugoslav Republic of Macedonia", "Togo", "Tokelau", "Tonga", "Trinidad and Tobago", "Tromelin Island", "Tunisia", "Turkey", "Turkmenistan", "Turks and Caicos Islands", "Tuvalu", "USA", "Uganda", "Ukraine", "United Arab Emirates", "United Kingdom", "Uruguay", "Uzbekistan", "Vanuatu", "Venezuela", "Viet Nam", "Virgin Islands", "Wake Island", "Wallis and Futuna", "West Bank", "Western Sahara", "Yemen", "Zambia", "Zimbabwe"}

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
    "WWTP" : {"organism" : "uncultured crAssphage", "wastewater_type" : "human waste", "sewage_type" : "municiple"},
    'fecal sample': {"organism" : "uncultured crAssphage"},
    'post-STP, refugee camp': {"organism" : "uncultured crAssphage", "wastewater_type" : "human waste", "sewage_type" : "municiple"},
    'pre-sewage treatment plant': {"organism" : "uncultured crAssphage", "wastewater_type" : "human waste", "sewage_type" : "municiple"},
    'Raw Sewage': {"organism" : "uncultured crAssphage", "wastewater_type" : "human waste", "sewage_type" : "municiple"}
}

columns = {
    "date" : "collection_date", "lat_lon" : "lat_lon", "name" : "sample_name",
    "address" : "adress", "altitude" : "altitude", "method" : "extraction_method", "sex" : "host_sex",
    "volunteer" : "host_subject_id", "locality" : "locality", "note" : "note", "contact" : "provider",
    "sampletype" : "sample_frequency", "description" : "sample_title", "source" : "source",
    "university" : "university"
}

def parse_file(filename):
    """
    PArse a tsv file
    :param filename: what is the name
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
            data[p[0]] = {}

            # deal with the country
            if "country" not in headers:
                sys.stderr.write(f"{bcolors.FAIL}FATAL:{bcolors.ENDC} No country found in {f}\n")
                sys.exit(-1)
            cidx = headers.index("country")
            if p[cidx] and p[cidx] not in countries:
                sys.stderr.write(f"{bcolors.FAIL}FATAL:{bcolors.ENDC}: |{p[cidx]}| is not a valid country\n")
                sys.exit(-1)
            if p[cidx]:
                data[p[0]]['geo_loc_name'] = p[cidx]


            # deal with the source
            src = p[sourceidx]
            if "SRA" == src:
                # Uncultivated Euryarchaeota archaeon UBA41 genome recovered from ERX556009
                bpidx = headers.index('BioProject')
                data[p[0]]['sample_name'] = f"Uncultured crAssphage amplicon recovered from {p[bpidx]}"
            else:
                # deal with envO terms
                if src not in env_broad_scale:
                    sys.stderr.write(f"{bcolors.FAIL}FATAL:{bcolors.ENDC}: |{src}| is not a valid source\n")
                    sys.exit(-1)
                data[p[0]]['env_broad_scale'] = env_broad_scale[src]
                data[p[0]]['env_local_scale'] = env_local_scale[src]
                data[p[0]]['env_medium'] = env_medium[src]

                # predefined terms
                for pd in predefined[src]:
                    data[p[0]][pd] = predefined[src][pd]

            # now the rest:
            for c in columns:
                if c not in headers:
                    sys.stderr.write(f"{bcolors.WARNING}WARN:{bcolors.ENDC} {c} not found in {f}\n")
                    continue
                idx = headers.index(c)
                if p[idx]:
                    data[p[0]][columns[c]] = p[idx]
            if 'description' not in data[p[0]]:
                data[p[0]]['description'] = data[p[0]]['sample_name']

    return data


def print_all(data):
    """
    Print out everything
    :param data:
    :return:
    """

    order = ["sample_name", "sample_title", "bioproject_accession", "organism", "collection_date", \
            "env_broad_scale", "env_local_scale", "env_medium", "geo_loc_name", "host", "lat_lon"]

    # figure out all the keys
    ak = set()
    for d in data:
        ak.update(data[d].keys())
    for k in ak:
        if k in order:
           ak.remove(k)
    otherkeys = sorted(ak)

    for d in data:
        sys.stdout.write(d)
        for o in order:
            if o not in data[d]:
                data[d][o] = ""
            sys.stdout.write(f"\t{data[d][o]}")
        for o in otherkeys:
            if o not in data[d]:
                data[d][o] = ""
            sys.stdout.write(f"\t{data[d][o]}")
        sys.stdout.write("\n")




__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse files for the NCBI BioSample Submission')
    parser.add_argument('-d', help='directory of csv files', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    data = {}
    for f in os.listdir(args.d):
        sys.stderr.write(f"{bcolors.OKGREEN}PARSING:{bcolors.ENDC} {f}\n")
        newdata = parse_file(os.path.join(args.d, f))
        data.update(newdata)

    print_all(data)
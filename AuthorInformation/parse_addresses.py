"""

"""

import os
import sys
import argparse
from author import Author, Address
import operator
from collections import OrderedDict

def parse_file(filename):
    """
    parse a file and create a set of authors
    :param filename: file to parse
    :return: set of author elements
    """

    authors = set()
    firstline = True # ignore the first line

    with open(filename, 'r') as f:
        for l in f:
            if firstline:
                firstline = False
                continue
            p = l.rstrip().split("\t")

            if len(p) < 15:
                sys.stderr.write("ERROR: Malformed: {}\t{}\n".format(len(p), p))
                continue

            auth = Author(p[2])

            try:
                if p[1]:
                    auth.orcid = p[1]
                if p[3]:
                    auth.lastname = p[3]
                    auth.lastnamelower = p[3].lower()
                if p[4]:
                    auth.firstname = p[4]
                    auth.firstnamelower = p[4].lower()
                if p[5]:
                    auth.middleinitial = p[5]
                if p[6]:
                    auth.email = p[6]
                if p[14]:
                    auth.order = int(p[14])
                if p[15]:
                    auth.contribution = p[15]

                primary = Address()
                if p[7]:
                    primary.department = p[7]
                if p[8]:
                    primary.institution = p[8]
                if p[9]:
                    primary.street = p[9]
                if p[10]:
                    primary.city = p[10]
                if p[11]:
                    primary.state = p[11]
                if p[12]:
                    primary.zip = p[12]
                if p[13]:
                    primary.country = p[13]

                auth.primaryaddress = primary

                secondary = Address()
                if len(p) > 17 and p[17]:
                    secondary.department = p[17]
                if len(p) > 18 and p[18]:
                    secondary.institution = p[18]
                if len(p) > 19 and p[19]:
                    secondary.street = p[19]
                if len(p) > 20 and p[20]:
                    secondary.city = p[20]
                if len(p) > 21 and p[21]:
                    secondary.state = p[21]
                if len(p) > 22 and p[22]:
                    secondary.zip = p[22]
                if len(p) > 23 and p[23]:
                    secondary.country = p[23]

                if secondary.is_valid():
                    auth.secondaryaddress = secondary
            except Exception as err:
                sys.stderr.write("Error parsing {}: {}\n".format(p[2], err))
                continue

            authors.add(auth)
    return authors


def test_validity(authors):
    """
    Test whether the author information is valid
    :param authors: the set of authors
    :return: nothing
    """

    abbs = set()

    for a in authors:
        if a.abbreviation.lower() in abbs:
            sys.stderr.write("FATAL: {} is a duplicate abbrevation\n".format(a.abbreviation))
        abbs.add(a.abbreviation.lower())
        if not a.is_valid():
            a.verbose = True
        print("{}\t{}\t{}".format(a.abbreviation, a.get_name(), a.is_valid()))

def check_spellings(authors):
    """
    Check for potential misspellings based on institutions, departments, and zip codes
    :param authors: the list of authors
    :return:
    """

    # find the set of different addresses
    addresses = {}
    allinst = {}
    allzip = {}
    for a in sorted(authors, key=operator.attrgetter('order', 'lastnamelower', 'firstnamelower')):
        pa = a.primaryaddress.get_address()
        if pa not in addresses:
            addresses[pa] = {}
            addresses[pa]['institution'] = a.primaryaddress.institution
            if a.primaryaddress.institution not in allinst:
                allinst[a.primaryaddress.institution] = 1
            else:
                allinst[a.primaryaddress.institution] += 1
            addresses[pa]['zip'] = a.primaryaddress.zip
            if a.primaryaddress.zip not in allzip:
                allzip[a.primaryaddress.zip] = 1
            else:
                allzip[a.primaryaddress.zip] += 1

        if a.secondaryaddress.is_valid():
            sa = a.secondaryaddress.get_address()
            if sa not in addresses:
                addresses[sa] = {}
                addresses[sa]['institution'] = a.secondaryaddress.institution
                if a.secondaryaddress.institution not in allinst:
                    allinst[a.secondaryaddress.institution] = 1
                else:
                    allinst[a.secondaryaddress.institution] += 1
                addresses[sa]['zip'] = a.secondaryaddress.zip
                if a.secondaryaddress.zip not in allzip:
                    allzip[a.secondaryaddress.zip] = 1
                else:
                    allzip[a.secondaryaddress.zip] += 1

    sys.stderr.write("Duplicates by institution\n")
    for i in allinst:
        if not i:
            continue
        if allinst[i] < 2:
            continue
        sys.stderr.write("\n")
        for a in addresses:
            if addresses[a]['institution'] == i:
                sys.stderr.write("{}\n".format(a))

    sys.stderr.write("\n\n\nDuplicates by zip\n")
    for z in allzip:
        if not z:
            continue
        if allzip[z] < 2:
            continue
        sys.stderr.write("\n")
        for a in addresses:
            if addresses[a]['zip'] == z:
                sys.stderr.write("{}\n".format(a))


def print_author_list(authors):
    """
    Print the list of all authors.
    :param authors: the set of authors
    :return:
    """

    # the list of addresses as we add them. This becomes the order
    addresses = []
    a: Author
    sys.stdout.write("<p>")
    for a in sorted(authors, key=operator.attrgetter('order', 'lastnamelower', 'firstnamelower')):
        output = a.get_name()
        pa = a.primaryaddress.get_address()
        if pa not in addresses:
            addresses.append(pa)
        addidx = addresses.index(pa) + 1
        if a.secondaryaddress.is_valid():
            sa = a.secondaryaddress.get_address()
            if sa not in addresses:
                addresses.append(sa)
            oidx = addresses.index(sa) + 1
            addidx = "{},{}".format(addidx, oidx)
        output += "<sup>{}</sup>, ".format(addidx)
        sys.stdout.write(output)

    sys.stdout.write("</p>\n\n\n")
    for i, j in enumerate(addresses):
        print("<sup>{}</sup>{}<br>".format(i+1,j))

def print_author_contributions(authors):
    """
    Print the author contribution list
    :param authors:
    :return:
    """

    contribs = OrderedDict()


    a: Author
    for a in sorted(authors, key=operator.attrgetter('order', 'lastnamelower', 'firstnamelower')):
        c = a.contribution
        thisc = ''.join([c[0].lower(), c[1:]]) # convert the first letter to lower case as it will be in a sentence
        if thisc not in contribs:
            contribs[thisc] = []
        contribs[thisc].append(a.abbreviation)

    sys.stdout.write("<p> &nbsp; </p><h1>Author Contributions</h1><p> &nbsp; </p>\n")
    for c in contribs:
        output = ", ".join(map(str, sorted(contribs[c])))
        output += " {}".format(c)
        sys.stdout.write("{}. ".format(output))
    sys.stdout.write("</p>\n\n\n")


def nature_form(authors):
    """
    Print a text list of authors to be cut and pasted for nature

    First Name
    Middle Initial
    Last Name
    email
    Institution
    City
    Country

    :param authors: the set of authors
    :return:
    """
    counter=0
    for a in sorted(authors, key=operator.attrgetter('order', 'lastnamelower', 'firstnamelower')):
        if None == a.middleinitial:
            a.middleinitial = ""
        print(counter)
        print("\n".join(map(str, [a.firstname, a.middleinitial, a.lastname, a.email,
                         a.primaryaddress.institution, a.primaryaddress.city,
                         a.primaryaddress.country])))
        print("\n\n")
        counter+=1

def science_list(authors):
    """
    Print a list of the authors for science. This is easy, first, last, email
    :param authors: the set of authors
    :return:
    """
    for a in sorted(authors, key=operator.attrgetter('order', 'lastnamelower', 'firstnamelower')):
        print("\t".join(map(str, [a.firstname, a.lastname, a.email])))


def email_list(authors):
    """
    Print a list suitable for emailing all colleagues.
    :param authors: the set of authors
    :return:
    """
    for a in sorted(authors, key=operator.attrgetter('order', 'lastnamelower', 'firstnamelower')):
        sys.stdout.write("{} {} <{}>, ".format(a.firstname, a.lastname, a.email))
    print()


def orcid_list(authors):
    """
    Print the author list in the correct order, but print their orcid
    :param authors: the set of authors
    :return:
    """

    for a in sorted(authors, key=operator.attrgetter('order', 'lastnamelower', 'firstnamelower')):
        orcid = a.orcid
        nm = a.get_name()
        if orcid:
            print(f"{orcid}\t{nm}")
        else:
            sys.stderr.write("No ORCID for {}\n".format(nm))
            print(nm)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Parse the author information from our google doc for the crAssphage paper")
    parser.add_argument('-f', help='Google doc of author information', required=True)
    parser.add_argument('-t', help="test validity of the author information", action='store_true')
    parser.add_argument('-d', help='check for duplicate entries', action='store_true')
    parser.add_argument('-o', help='print the author list as ORCids in the correct order', action='store_true')
    parser.add_argument('-n', help='print the author list suitable for cutting and pasting to nature', action='store_true')
    parser.add_argument('-s', help='print the author list to add to the science bulk upload', action='store_true')
    parser.add_argument('-e', help='print the author list to use sending emails', action='store_true')
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    authors = parse_file(args.f)

    if args.d:
        check_spellings(authors)
        sys.exit(0)

    if args.t:
        test_validity(authors)
        sys.exit(0)

    if args.n:
        nature_form(authors)
        sys.exit(0)

    if args.s:
        science_list(authors)
        sys.exit(0)

    if args.e:
        email_list(authors)
        sys.exit(0)

    if args.o:
        orcid_list(authors)
        sys.exit(0)

    print_author_list(authors)
    print_author_contributions(authors)




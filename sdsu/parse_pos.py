from bs4 import BeautifulSoup
import re
import sys


if len(sys.argv) < 2:
    sys.exit("{} <file to parse> <verbose output>".format(sys.argv[0]))

verbose = False
if len(sys.argv) > 2:
    verbose = True

#soup = BeautifulSoup(open(sys.argv[1], 'r'), 'lxml')
soup = BeautifulSoup(open(sys.argv[1], 'r'), 'html.parser')

wanted = {"Student ID", "First Name", "Middle Name", "Last Name", "Address", "City", "State", "Zip Code",
                            "Country", "Phone1", "Phone2", "E-mail", "SIMS Major Code", "HEGIS Major Code", "Major", "Degree", 
                            "Student Standing", "Term"}


for s in soup.body.find_all(text=re.compile('Student\s*Summary')):
    t = s.find_parent('table')
    if verbose:
        sys.stderr.write("Found a table for Student Summary\n")
    for row in t.find_all("tr")[1:]:
        r = [re.sub('\s+', ' ', cell.get_text(strip=True)) for cell in row.find_all("td")]
        if r[0] in wanted:
            print("{}\t{}".format(r[0], r[1]))
        if len(r) > 2:
            if r[2] in wanted:
                print("{}\t{}".format(r[2], r[3]))

print("CLASSES")

#if verbose:
#    print(soup.body)

# this method doesn't work now. For some reason the SDSU Course table is just the header.
if False:
    for s in soup.find_all(text=re.compile('SDSU\s*Course')):
        t = s.find_parent('table')
        if verbose:
            print("\nFound a table for Student Course:\n>>>\n")
            print(str(t) + "\n<<<<\n")
        for row in t.find_all("tr")[1:]:
            print("\t".join([re.sub('\s+', ' ', cell.get_text(strip=True)) for cell in row.find_all("td")]))


for s in soup.find_all("tr"):
    cells = s.find_all("td")
    if len(cells) == 9:
        x = "\t".join([re.sub('\s+', ' ', c.get_text(strip=True)) for c in cells])
        if "wrapperMain" in x:
            continue
        if "Department" in x:
            continue
        print(x)

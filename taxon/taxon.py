import gzip

# the location of the taxonomy files
defaultdir = '/home/db/taxonomy/current/'

'''
From nodes.dmp
        tax_id                                  -- node id in GenBank taxonomy database
        parent tax_id                           -- parent node id in GenBank taxonomy database
        rank                                    -- rank of this node (superkingdom, kingdom, ...)
        embl code                               -- locus-name prefix; not unique
        division id                             -- see division.dmp file
        inherited div flag  (1 or 0)            -- 1 if node inherits division from parent
        genetic code id                         -- see gencode.dmp file
        inherited GC  flag  (1 or 0)            -- 1 if node inherits genetic code from parent
        mitochondrial genetic code id           -- see gencode.dmp file
        inherited MGC flag  (1 or 0)            -- 1 if node inherits mitochondrial gencode from parent
        GenBank hidden flag (1 or 0)            -- 1 if name is suppressed in GenBank entry lineage
        hidden subtree root flag (1 or 0)       -- 1 if this subtree has no sequence data yet
        comments                                -- free-text comments and citations
'''


class taxonNode:
    def __init__(self, t = None, p = None, r = None, e = None, d = None, i = None, gc = None, igc = False, mgc = None, imgc = False,  gh = False, hs = False, c=None, *others):
        self.parent = p
        self.taxid = t
        self.rank = r
        self.embl = e
        self.division = d
        self.inherited = i
        self.geneticCode = gc
        self.inheritedGC = igc
        self.mitochondrialGeneticCode = mgc
        self.inheritedMitochondrialGeneticCode = imgc
        self.GenBankHidden = gh
        self.hiddenSubtree = hs
        self.comments = c
        if len(others) > 0:
            print "WARNING: ", p, " :: ", others

'''
Taxonomy names file (names.dmp):
        tax_id                                  -- the id of node associated with this name
        name_txt                                -- name itself
        unique name                             -- the unique variant of this name if name not unique
        name class                              -- (synonym, common name, ...)
'''

class taxonName:
     def __init__(self, t = None, n = None, u = None, nc = None):
        self.taxid = t
        self.name = n
        self.unique = u
        self.nameClass = nc
        
'''
Divisions file (division.dmp):
        division id                             -- taxonomy database division id
        division cde                            -- GenBank division code (three characters)
        division name                           -- e.g. BCT, PLN, VRT, MAM, PRI...
        comments
'''
class taxonDivision:
     def __init__(self, i = None, c = None, n = None, co = None):
        self.divid = i
        self.name = n
        self.code = c
        self.comments = co
        




def readTaxa():
    '''Read the taxonomy tree. An alias for readNodes()'''
    return readNodes()


def readNodes():
    '''Read the node information from the default location'''
    taxa = {}
    fin = open(defaultdir+'nodes.dmp', 'r')
    for line in fin:
        line = line.rstrip('\t|\n')
        cols = line.split('\t|\t')
        t=taxonNode(*cols)
        taxa[cols[0]]=t
    fin.close()
    return taxa


def extendedNames():
    '''
    Extended names returns "genbank synonym" and "synonym" as well as 
    "scientific name" and "blast name". Because we are reading more 
    names it is slower and consumes more memory
    '''
    names = {}
    blastname={}
    genbankname={}
    synonym={}
    fin = open(defaultdir+'names.dmp', 'r')
    for line in fin:
        line = line.rstrip('\t|\n')
        cols = line.split('\t|\t')
        t=taxonName(*cols)
        if "scientific name" in cols[3]:
            names[cols[0]]=t
        elif "blast name" in cols[3]:
            blastname[cols[0]]=t
        elif "genbank synonym" in cols[3]:
            genbankname[cols[0]]=t
        elif "synonym" in cols[3]:
            synonym[cols[0]]=t

    fin.close()
    return names, blastname, genbankname, synonym



def readNames():
    '''Read the name information from the default location'''
    names = {}
    blastname={}
    fin = open(defaultdir+'names.dmp', 'r')
    for line in fin:
        line = line.rstrip('\t|\n')
        cols = line.split('\t|\t')
        t=taxonName(*cols)
        if "scientific name" in cols[3]:
            names[cols[0]]=t
        if "blast name" in cols[3]:
            blastname[cols[0]]=t
    fin.close()
    return names, blastname

def readDivisions():
    '''Read the divisions.dmp file'''
    divs = {}
    fin = open(defaultdir+'division.dmp', 'r')
    for line in fin:
        line = line.rstrip('\t|\n')
        cols = line.split('\t|\t')
        t=taxonDivision(*cols)
        divs[cols[0]]=t
    fin.close()
    return divs


def readGiTaxId(type='nucl'):
    '''
    Read gi_taxid.dmp. You can specify the type of database that you
    want to parse, default is nucl (nucleotide), can also accept prot 
    (protein).
    
    Returns a hash of gi and taxid
    '''
    if type != 'nucl' and type != 'prot':
        sys.stderr.write("Type must be either nucl or prot, not " + type + "\n")
        sys.exit(-1)
    fileIn = defaultdir + "/gi_taxid_" + type + ".dmp.gz"
    taxid={}
    fin = gzip.open(fileIn, 'r')
    for line in fin:
        line = line.strip()
        parts=line.split("\t")
        taxid[parts[0]]=parts[1]
    fin.close()
    return taxid


def readTaxIdGi(type='nucl'):
    '''
    Read gi_taxid.dmp. You can specify the type of database that you
    want to parse, default is nucl (nucleotide), can also accept prot 
    (protein).
    
    NOTE: This method returns taxid -> gi not the other way around. This
    may be a one -> many mapping (as a single taxid maps to more than
    one gi), and so we return a list of gi's for each taxid.

    Returns a hash of taxid and gi
    '''
    if type != 'nucl' and type != 'prot':
        sys.stderr.write("Type must be either nucl or prot, not " + type + "\n")
        sys.exit(-1)
    fileIn = defaultdir + "/gi_taxid_" + type + ".dmp.gz"
    taxid={}
    fin = gzip.open(fileIn, 'r')
    for line in fin:
        line = line.strip()
        parts=line.split("\t")
        if parts[1] not in taxid:
            taxid[parts[1]]=[]
        taxid[parts[1]].append(parts[0])
    fin.close()
    return taxid



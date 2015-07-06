import sys
from servers.SAP import SAPserver

server=SAPserver()
genomes=server.all_genomes({'-complete':1, '-prokaryotic':1})
print("There are ", len(genomes), " genomes")

#sys.exit(0)
for genomeID in genomes:
    sys.stderr.write("Genome: " + genomes[genomeID] + " (" + str(genomeID) + ")\n")
    prots = server.all_proteins({'-id':genomeID})
    print("There are ", len(prots), " prots")
    dna   = server.ids_to_sequences({'-ids':list(prots.keys())})
    print("There are ", len(dna), " sequences")
    funcs = server.ids_to_functions({'-ids':list(prots.keys())})
    print("There are ", len(funcs), " functions")


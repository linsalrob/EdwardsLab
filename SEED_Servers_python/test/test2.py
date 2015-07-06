
import sys
from servers.SAP import SAPserver

server=SAPserver()


genomeID = '83333.1'
sys.stderr.write("Genome: " + str(genomeID) + "\n")
prots = server.all_proteins( {"-id" : genomeID} )
print("protein length " + str(len(prots)))

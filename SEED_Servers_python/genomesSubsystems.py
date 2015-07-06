
import sys
from servers.SAP import SAPserver

server=SAPserver()
genomes = server.all_genomes({'-complete':1, '-prokaryotic':1})
allss = server.all_subsystems(None)
ss = list(allss.keys())
ss.sort()
subsystems = server.genomes_to_subsystems( {'-ids':list(genomes.keys()), '-all':1}  )

gHasSS={}
for g in subsystems:
    gHasSS[g]={}
    for tpl in subsystems[g]:
        if tpl[1] == '0':
            continue
        gHasSS[g][tpl[0]]=tpl[1]

print("", "", *ss, sep='\t')
for g in genomes:
    print(genomes[g], g, sep='\t', end="")
    for s in ss:
        if s in gHasSS[g]:
            print("\t", gHasSS[g][s], end="")
        else:
            print("\t", end="")
    print()





from Bio import SeqIO

seq = SeqIO.read('sequence.gb', 'genbank')
print(seq.id)


with open('features.tsv', 'w') as out:
    for feature in seq.features:
        if 'locus_tag' in feature.qualifiers:
            lt = feature.qualifiers['locus_tag'][0]
        if 'product' in feature.qualifiers:
            out.write(lt + "\t" + feature.qualifiers['product'][0] + "\n")
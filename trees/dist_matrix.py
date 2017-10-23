import sys
import itertools
from ete3 import Tree

try:
    t = Tree()
    t.populate(int(sys.argv[1]), random_branches=True)
except ValueError:
    print >>sys.stderr, 'loading', sys.argv[1]
    t = Tree(sys.argv[1])

lineages = {}
for tip in t:
    lin = []
    n = tip
    while n.up:
        lin.append(n)
        n = n.up
    lineages[tip.name] = set(lin)

matrix = {}
def get_dist(a, b):
    if a == b:
        return 0.0
    try:
        return matrix[(a, b)]
    except KeyError:
        return matrix[(b, a)]

for tip_a, tip_b in itertools.permutations(lineages.keys(), 2):
    d = sum([n.dist for n in lineages[tip_a] ^ lineages[tip_b]])
    matrix[(tip_a, tip_b)] = d
    #if len(matrix) % 10000 == 0:
    #    print >>sys.stderr, len(matrix)

leaves = t.get_leaf_names()
print '\t'.join(['#names'] + leaves)
for tip_a in leaves:
    row = [tip_a]
    for tip_b in leaves:
        row.append(get_dist(tip_a, tip_b))
    print '\t'.join(map(str, row))


# test

import random
s = random.sample(matrix.keys(), 1000)
for a,b in s:
    d0 = get_dist(a, b)
    d1 = t.get_distance(a, b)
    if round(d0, 8) != round(d1, 8):
        print >>sys.stderr, a, b, d0, d1
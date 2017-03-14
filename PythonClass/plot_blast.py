import matplotlib.pyplot as plt
import re
import sys

longest_sequence = 1601

minlen = 45

hits = {}

with open('/home/redwards/Desktop/16S_all_padded.blastn', 'r') as f:
    for l in f:
        p = l.strip().split('\t')
        for i in [2,10,11]:
            p[i] = float(p[i])
        for i in range(3,10):
            p[i] = int(p[i])

        if p[8] > p[9]:
            (p[8], p[9]) = (p[9], p[8])

        m = re.match('fig\|(\d+\.\d+)', p[1])
        if not m:
            sys.stderr.write("Can't parse " + p[1] + "\n")
            continue
        sid = m.group(1)

        if sid not in hits:
            hits[sid]=[0] * longest_sequence

        for i in range(p[8], p[9]+1):
            hits[sid][i] = hits[sid][i] + 1


fig = plt.figure()
ax = fig.add_subplot(111)


x=range(longest_sequence)

for sid in hits:
    ax.plot(x, hits[sid], label=sid)

ax.legend()
fig.set_facecolor('white')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

plt.show()

















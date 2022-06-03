import sys

assert len(sys.argv) == 2
genome = sys.argv[1]
print('genome:', genome)

suffixes = [genome[i:] for i in range(len(genome))]
order = {suf: i for i, suf in enumerate(suffixes)}

suffixes = sorted(suffixes)
sa = [order[suf] for suf in suffixes]

for e in sa:
    print(e, end=' ')
print()
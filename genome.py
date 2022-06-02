genome = "ACGTACACACCCGCTACCGACCGTC$"

suffixes = [genome[i:] for i in range(len(genome))]
order = {suf: i for i, suf in enumerate(suffixes)}

suffixes = sorted(suffixes)
sa = [order[suf] for suf in suffixes]

print(sa)
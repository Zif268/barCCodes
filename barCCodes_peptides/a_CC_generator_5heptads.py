##############
### IMPORT ###
##############

import os
import sys
import pandas as pd
import time
from itertools import permutations

#################
### CONSTANTS ###
#################

# heptad length
h_l: int = 5

# export
origin = sys.argv[1]
dirname = os.path.join(origin, "Sequences")
os.makedirs(dirname, exist_ok=True)
out_file = os.path.join(dirname, "CC_list_all_" + str(h_l) + "heptads_AQS.csv")

# heptad options
heptads = ["EI", "EN", "KI", "KN"]

bcf = ["A", "A", "A", "Q", "Q", "Q", "S", "S", "S"]
bcf_comb = list(sorted(set(permutations(bcf, 3))))

# empty list for generation of heptad combinations
pep = []

###################
### DEFINITIONS ###
###################


def generate_cc(p):
    """
    generates all possible combinations of heptads for the defined peptide length.
    :param p: empty list
    :return: list of lists with all possible heptad combinations
    (example: [["EN", "KI", "KI", "EI", "EN"], ...]
    """
    global pep
    if len(p) == h_l:
        pep = pep + [p]
        return
    else:
        for h in heptads:
            generate_cc(p + [h])


def heptad_seq(h):
    """
    generates all possible amino acid sequences for heptads with permutations at b, c and f positions
    :param h: string representing a heptad (example: "KN")
    :return: list of 7-residue sequences (strings) for all permutations (variants) of the given heptad
    (example: ["KNAALKA", "KNAALKQ", "KNAQLKQ", ..., "KNSSLKS"])
    """
    h_lst = []
    for pos in bcf_comb:
        h_seq = h + pos[0] + pos[1] + "L" + h[0] + pos[2]
        h_lst.append(h_seq)
    return h_lst


def peptide_seq(h):
    """
    generates full peptide sequences for all permutations of a heptad combination
    :param h: list of heptads (example: ["EN", "KI", "KI", "EI", "EN"])
    :return: dictionary with permutations (variants) as keys and full peptide sequences as values
    (example: {"AAA": SPEDENAALEAKIAALKAKIAALKAEIAALKAENAALEAG, "AAQ": SPEDENAALEQKIAALKQKIAALKQEIAALKQENAALEQG, ...}
    """
    sequences = dict()
    full_seqs = []
    variant_dict = dict()
    for j, hep in enumerate(h):
        sequences[j + 1] = heptad_seq(hep)
    for m in list(zip(sequences.get(1), sequences.get(2), sequences.get(3), sequences.get(4), sequences.get(5))):
        full_seqs.append("SPED" + str("".join(m)) + "G")
    for se in full_seqs:
        variant_dict[str(se[6]) + str(se[7]) + str(se[10])] = se
    return variant_dict


####################
### GENERATE CCs ###
####################


start_time = time.time()

# generate a list of all possible heptad combinations
generate_cc([])

# generate list of peptide names
pep_names = []
for i in range(1, len(pep) + 1):
    pep_names.append(i)

# dictionary of heptads with assigned names
# {1: ["EI", "EI", "EI", "EI", "EI"], 2: ["EI", "EI", "EI", "EI", "EN"], ...}
pep_d = dict(zip(pep_names, pep))

# generate a full dictionary with all permutations (variants) for all generated peptides
# {("name, variant"): [["EN", "KN", "KI", "EI", "EN"], "SPEDENAALEAKNAALKAKIAALKAEIAALEAENAALEAG"], ...}
full_dict = dict()
for k, v in pep_d.items():
    d = peptide_seq(v)
    for k1, v1 in d.items():
        full_dict[k, k1] = [v, v1]


##################
### WRITE FILE ###
##################


# create empty row for each CC and fill it
n = 0
rows = []
for k, v in full_dict.items():
    n = n + 1
    r = {"Peptide_name": k,
         "Number": list(k)[0],
         "Variant": list(k)[1],
         "Heptad1": v[0][0],
         "Heptad2": v[0][1],
         "Heptad3": v[0][2],
         "Heptad4": v[0][3],
         "Heptad5": v[0][4],
         "Sequence": v[1]
         }
    rows.append(r)

# write full list
table = pd.DataFrame(rows, index=None, columns=["Peptide_name",
                                                "Number",
                                                "Variant",
                                                "Heptad1",
                                                "Heptad2",
                                                "Heptad3",
                                                "Heptad4",
                                                "Heptad5",
                                                "Sequence"])
table.to_csv(out_file, header=True, index=None)
print(str(len(table)) + " peptides written.")


# generate separate csv files for each set
sets = set()
for i in list(full_dict.keys()):
    sets.add(i[1])
s = list(sorted(sets))
for i in s:
    seqs = dict()
    for k, v in full_dict.items():
        if k[1] == i:
            seqs[k[0]] = v
    out_f = os.path.join(dirname, "CC_list_AQS_" + str(h_l) + str(i) + ".csv")
    n = 0
    rows = []
    for x, y in seqs.items():
        n = n + 1
        r = {"Peptide_name": x,
             "Heptad1": y[0][0],
             "Heptad2": y[0][1],
             "Heptad3": y[0][2],
             "Heptad4": y[0][3],
             "Heptad5": y[0][4],
             "Sequence": y[1]
             }
        rows.append(r)
    t = pd.DataFrame(rows, index=None, columns=["Peptide_name",
                                                "Heptad1",
                                                "Heptad2",
                                                "Heptad3",
                                                "Heptad4",
                                                "Heptad5",
                                                "Sequence"])
    t.to_csv(out_f, header=True, index=None)


print("--- %s seconds ---" % (time.time() - start_time))

##############
### Import ###
##############
import os
import sys
import pandas as pd
from itertools import permutations
import time
import math

#################
### Constants ###
#################

h_l = 4

# generate a csv file with interaction scores for all sets
origin = sys.argv[1]
indir = os.path.join(origin, "Helicities")
outdir = os.path.join(origin, "Interaction_scores")
os.makedirs(outdir, exist_ok=True)
out_file_all = os.path.join(outdir, "CC_4heptads_all_scores-v1.0.csv")
df_all = pd.DataFrame()
# print(df_all)

# Compute all permutations
bcf = ["A", "A", "A", "Q", "Q", "Q", "S", "S", "S"]
bcf_comb = list(sorted(set(permutations(bcf, 3))))


# coupling energies for each heptad combination in kcal/mol
coup_en = {("EI", "KI"): -2.7,
           ("KI", "EI"): -2.7,
           ("EN", "KN"): -2.3,
           ("KN", "EN"): -2.3,
           ("EI", "EN"): 5.4,
           ("EN", "EI"): 5.4,
           ("KI", "KN"): 5.2,
           ("KN", "KI"): 5.2,
           ("KN", "EI"): 3.1,
           ("EI", "KN"): 3.1,
           ("EN", "KI"): 2.5,
           ("KI", "EN"): 2.5,
           ("EI", "EI"): -0.1,
           ("KI", "KI"): -0.3,
           ("EN", "EN"): 0.3,
           ("KN", "KN"): 0.1,
           }

# gas constant R for calculation of Kd in kcal/(K*mol)
r_const = 0.00198720425864083

# temperature in K (310.15 = 37°C)
temperature = 310.15

# correct pairings of heptads
hep_pairs = [("EN", "KN"), ("KN", "EN"), ("KI", "EI"), ("EI", "KI")]

# types of interactions
int_type = {"pair": [4, 4],
            "shift_0:mismatch_1": [4, 3],
            "shift_0:mismatch_2": [4, 2],
            "shift_0:mismatch_3": [4, 1],
            "shift_0:mismatch_4": [4, 0],
            "shift_1:mismatch_0": [3, 3],
            "shift_1:mismatch_1": [3, 2],
            "shift_1:mismatch_2": [3, 1],
            "shift_1:mismatch_3": [3, 0],
            "shift_2:mismatch_0": [2, 2],
            "shift_2:mismatch_1": [2, 1],
            "shift_2:mismatch_2": [2, 0]}

############################
### Function definitions ###
############################


def conf(p1, p2):
    """
    finds all possible pairing conformations for two 4heptad peptides up to two interacting heptads
    :param p1: list of heptads representing a peptide (example: ["EI", "EI", "EI", "EI"])
    :param p2: list of heptads representing a peptide (example: ["EI", "EI", "EI", "EI"])
    :return: tuple of tuples with possible heptad pairings
    (example: ((("EI", "KI"), ("KI, "EN), ("EN", "KN")), (("EI, "EI"), ("KN", "KI"))))
    """
    z0 = tuple(zip(p1, p2))
    z1l = tuple(zip(p1, p2[1:]))
    z1r = tuple(zip(p1[1:], p2))
    z2l = tuple(zip(p1, p2[2:]))
    z2r = tuple(zip(p1[2:], p2))
    return z0, z1l, z1r, z2l, z2r


def dg_calc(t):
    """
   calculates the delta G for a given conformation of two paired peptides
   :param t: tuple with heptad pairings in the form (("EI", "KI"), ("KI, "EN), ("EN", "KN"), ("EI", "KI"))
   :return: sum of coupling energies (delta G)
   """
    dg = math.fsum(coup_en.get(h) for h in t)
    return round(dg, 1)


def calc_kd(dg):
    """
    calculates Kd from given delta G
    :param dg: delta G in kcal/mol (number)
    :return: calculated Kd (number)
    """
    return math.exp(dg / (r_const * temperature))


def pairing(p1, p2):
    """
    calculates coupling energies for all possible conformations of two peptides
    :param p1: list of heptads in the form ["EI", "EI", "EI", "EI"]
    :param p2: list of heptads in the form ["EI", "EI", "EI", "EI"]
    :return: two dictionaries with conformation tuples as keys and coupling energies (dg) and as values
    """
    conformations_dg = dict()
    for z in conf(p1, p2):
        conformations_dg[z] = dg_calc(z)
    return conformations_dg


def best_conf(cfs):
    """
    finds the coupling energy of the strongest conformation
    :param cfs: dictionary with conformation tuples as keys and coupling energies as values
    :return: lowest coupling energy, list of conformations with the lowest coupling energies
    """
    min_dg = min(cfs.values())
    min_conf = [cf for cf, deltag in cfs.items() if deltag == min_dg]
    return min_dg, min_conf


def int_types(min_c):
    """
    finds the interaction types of the best conformations
    :param min_c: list of conformations with lowest coupling energies
    :return: list of interaction types for conformations with lowest coupling energies
    """
    types = []
    for c in min_c:
        no_zip = len(c)
        no_match = 0
        for z in c:
            if z in hep_pairs:
                no_match = no_match + 1
        i_type = list(int_type.keys())[list(int_type.values()).index([no_zip, no_match])]
        types.append(i_type)
    return types


def calc_pf(hel1, hel2, diss):
    """
    calculates interaction score
    :param hel1: helicity of peptide 1 (number)
    :param hel2: helicity of peptide 2 (number)
    :param diss: Kd (number)
    :return: interaction score (number)
    """
    return (hel1 * hel2 * 0.000001) / diss


def is_pair(t):
    """
    determines if two peptides are a predicted pred_pair
    :param t: tuple of tuples with heptad combinations ([1] element in the output of best_conf function)
    :return: True if predicted pred_pair, False if not pred_pair
    """
    pred_pair = False
    tups = t[0]
    if len(tups) == h_l:
        pair_h = 0
        for h in tups:
            if h in hep_pairs:
                pair_h = pair_h + 1
        if pair_h == h_l:
            pred_pair = True
    return pred_pair


###############
### Process ###
###############

for comb in bcf_comb:
    ccset = "".join(comb)

    # import list of CCs with helicity data
    fname = os.path.join(indir, "CC_list_4heptads_" + str(ccset) + "_hel_filtered.csv")
    df = pd.read_csv(fname, header=0, index_col=0)

    peptides = dict()
    for i, row in df.iterrows():
        p = int(i.split(",")[0][1:])
        peptides[p] = [[row["Heptad1"], row["Heptad2"], row["Heptad3"], row["Heptad4"]],
                       row["Sequence"], row["Helicity"]]

    # generate a csv file
    out_file = os.path.join(outdir, "CC_4heptads" + str(ccset) + "_scores-v1.0.csv.gz")

    ########################
    ### Calculate scores ###
    ########################

    start_time = time.time()

    # dictionary with tuples of pairs as keys and kds and pairing factors as values
    no_pairs = 0
    kds = dict()
    for k1, v1 in peptides.items():
        for k2, v2 in peptides.items():
            best_con = best_conf(pairing(v1[0], v2[0]))[1]
            typ = int_types(best_con)
            pair = is_pair(best_con)
            if pair:
                no_pairs = no_pairs + 1
            kd = calc_kd(best_conf(pairing(v1[0], v2[0]))[0])
            int_score = calc_pf(v1[2], v2[2], kd)
            kds[tuple([k1, k2])] = [k1, k2, kd, pair, best_con, typ, int_score]

    df_kds = pd.DataFrame.from_dict(kds,
                                    orient='index',
                                    columns=['Peptide1', 'Peptide2', 'Kd', 'Pair', 'Best_conformation',
                                             'Interaction_type', 'Interaction_score'])
    cname = "IS_" + str(ccset)
    df_kds = df_kds.rename(columns={'Interaction_score': cname})

    if df_all.empty:
        df_all = df_kds
    else:
        df_all[cname] = df_kds[cname]

    # print(df_all)
    print(str(no_pairs) + " pairs found.")

    # create empty row for each combination and fill it
    n = 0
    rows = []
    for k, v in kds.items():
        n = n + 1
        r = {"Peptide1": k[0],
             "Peptide2": k[1],
             "Kd": v[2],
             "Best_conformation": v[4],
             "Interaction_type": v[5],
             "Interaction_score": v[6],
             "Pair": v[3],
             "P1_seq": peptides[k[0]][1],
             "P2_seq": peptides[k[1]][1]
             }
        rows.append(r)

    # write to disk
    table = pd.DataFrame(rows, index=None, columns=[
                                                    "Peptide1",
                                                    "Peptide2",
                                                    "Kd",
                                                    "Best_conformation",
                                                    "Interaction_type",
                                                    "Interaction_score",
                                                    "Pair",
                                                    "P1_seq",
                                                    "P2_seq"])
    table.to_csv(out_file, header=True, index=None, compression="gzip")
    print("Written %d rows to %s" % (len(table), out_file))

    print(str(len(table)) + " combinations written.")

    print("--- %s seconds ---" % (time.time() - start_time))

df_all.to_csv(out_file_all, header=True, compression="gzip")
print("Written %d rows to %s" % (len(df_all), out_file_all))

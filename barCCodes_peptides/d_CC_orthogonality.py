##############
### Import ###
##############
import sys

import pandas as pd
from itertools import permutations
import time
import os
import statistics
# import numpy as np

################
### Cut-offs ###
################


# cutoffs for orthogonal interaction
cutoff = [0.1, 0.2, 0.3, 0.4, 0.5]

# score for strong interaction
strong = 10

##############
### Inputs ###
##############
origin = sys.argv[1]
num = int(sys.argv[2])
model = sys.argv[3]

# Check heptad number
supported_heptads = (4, 5)
if num not in supported_heptads:
    raise NotImplementedError("Heptad number must be in %s" % ", ".join(list(map(str, supported_heptads))))

# Check supported model
supported_models = ("1.0", "2.1")
if model not in supported_models:
    raise NotImplementedError("Model must be in %s" % ", ".join(supported_models))


#############
### Files ###
#############

# table with screen data
rows_screen_all = []
seq_dir = os.path.join(origin, "Sequences")
hel_dir = os.path.join(origin, "Helicities")
int_dir = os.path.join(origin, "Interaction_scores")
dir_all = os.path.join(origin, "Orth_screen_all")
os.makedirs(dir_all, exist_ok=True)
out_screen_all = os.path.join(dir_all, "CC_%sheptads_orth-screen_all-v%s.csv" % (num, model))

bcf = ["A", "A", "A", "Q", "Q", "Q", "S", "S", "S"]
bcf_comb = list(sorted(set(permutations(bcf, 3))))


############################
### FUNCTION DEFINITIONS ###
############################

def orthogonal(pair1, pair2):
    """
    determines if two pairs are orthogonal
    :param pair1: peptide pair (tuple)
    :param pair2: peptide pair (tuple)
    :return: True if pairs are mutually orthogonal, False if not orthogonal
    """
    orth = True
    for pept1 in pair1:
        for pept2 in pair2:
            if tuple(sorted([pept1, pept2])) in non_pairs.keys():
                if non_pairs.get(tuple(sorted([pept1, pept2])))[1] >= orth_cutoff:
                    orth = False
    if orth:
        return True
    else:
        return False

###############
### Process ###
###############


for comb in bcf_comb:
    cc_set = "".join(comb)

    # calculate median helical propensity of the full peptide set
    fname = os.path.join(hel_dir, "CC_list_%dheptads_" % num + str(cc_set) + "_hel_filtered.csv")
    df_hp = pd.read_csv(fname, header=0)

    hps = []
    for i, row in df_hp.iterrows():
        hps.append(row["Helicity"])
    med_hp = statistics.median(hps)

    # import interaction scores
    fname = os.path.join(int_dir, "CC_%dheptads" % num + str(cc_set) + "_scores-v%s.csv.gz" % model)
    df = pd.read_csv(fname, header=0)

    peptides = dict()
    for i, row in df.iterrows():
        peptides[i] = [row["Peptide1"], row["Peptide2"], row["Kd"], row["Interaction_score"], row["Pair"]]

    # make dicts with pairs and non-pairs
    pairs = dict()
    non_pairs = dict()
    for k, v in peptides.items():
        if v[4]:
            pairs[tuple(sorted([v[0], v[1]]))] = [v[2], v[3]]
        if not v[4]:
            non_pairs[tuple(sorted([v[0], v[1]]))] = [v[2], v[3]]

    # import peptide sequences
    fname = os.path.join(seq_dir, "CC_list_AQS_%d" % num + str(cc_set) + ".csv")
    df2 = pd.read_csv(fname, header=0)
    seqs = dict()
    for i, row in df2.iterrows():
        seqs[row["Peptide_name"]] = row["Sequence"]

    # directory for export
    out_screen = os.path.join(dir_all, "CC_%dheptads" % num + str(cc_set) + "_orth-screen-v%s.csv" % model)


###########
### RUN ###
###########

    start_time = time.time()

    # table with screen data for individual sets
    rows_screen = []

    print(str(cc_set))

    # predict orthogonality and make a set of all possible sets of orthogonal pairs
    for orth_cutoff in cutoff:
        print(orth_cutoff)
        # find number of off-target interactions
        no_off = 0
        off_t = dict()
        for k, v in non_pairs.items():
            if non_pairs.get(k)[1] >= orth_cutoff:
                off_t[k] = v[1]
                no_off = no_off + 1
        print(str(no_off) + " off-target interactions found.")

        ort_sets = set()
        no_sets = 0
        hom = set()
        weak = 0

        ### GENERATE ALL POSSIBLE ORTHOGONAL SETS

        for k, v in pairs.items():

            ort_set = set()
            for pep in k:
                if non_pairs.get(tuple([pep, pep]))[1] >= orth_cutoff:
                    hom.add(pep)
                    break
                else:
                    if v[1] >= strong:
                        ort_set.add(k)
                    else:
                        weak = weak + 1

            if ort_set:
                no_sets = no_sets + 1

                for k1, v1 in pairs.items():

                    ort = True
                    for pep1 in k1:
                        if non_pairs.get(tuple([pep1, pep1]))[1] >= orth_cutoff:
                            ort = False
                            break
                        else:
                            if v1[1] <= strong:
                                ort = False
                                break
                            else:
                                for pep_pair in sorted(list(ort_set)):
                                    if not orthogonal(k1, pep_pair):
                                        ort = False
                                        break

                                if ort:
                                    ort_set.add(k1)

            pair_list = list()
            for key in list(ort_set):
                pair_list.append(key)
            ort_sets.add(tuple(sorted(pair_list)))

        print(str(len(hom)) + " peptides can form homodimers.")
        print(str(weak / 2) + " weak pairs found.")
        print("Generated " + str(no_sets) + " orthogonal peptide sets.")
        print("Sorted " + str(len(list(ort_sets))) + " different orthogonal peptide sets.")

        ### FIND LARGEST SETS
        max_length = max([len(i) for i in list(ort_sets)])
        largest_set = []
        for s in list(ort_sets):
            if len(s) == max_length:
                largest_set.append(s)

        print("The largest set contains " + str(max_length) + " orthogonal peptide pairs. " + str(len(largest_set)) +
              " sets of length " + str(max_length) + " found.")

        r_screen = {"Cutoff": "{:.2f}".format(orth_cutoff), "Off-target": no_off, "Homodimers": len(hom),
                    "Weak_pairs": weak / 2, "No_sets": len(list(ort_sets)),
                    "Largest_set": max_length, "Median_HP": med_hp}
        rows_screen.append(r_screen)

        r_screen_all = {"Set": cc_set, "Cutoff": "{:.2f}".format(orth_cutoff),
                        "Off-target": no_off, "Homodimers": len(hom), "Weak_pairs": weak / 2,
                        "No_sets": len(list(ort_sets)), "Largest_set": max_length, "Median_HP": med_hp}
        rows_screen_all.append(r_screen_all)

        ### SCORE LARGEST SETS

        # find best set by ratio of highest and lowest interaction scores and/or pair strength
        best_set = []

        if len(largest_set) > 1:
            setscore = dict()
            for s in largest_set:
                sscore = dict()
                for p in s:
                    sscore[p] = pairs.get(p)[1]
                maxp = max(zip(sscore.values(), sscore.keys()))[0]
                minp = min(zip(sscore.values(), sscore.keys()))[0]
                tot_str = sum(sscore.values())
                setscore[s] = [maxp / minp, tot_str]

            # find minimal ratio
            rats = []
            for i in list(setscore.values()):
                rats.append(i[0])
            min_rat = min(rats)

            # filter set by minimal ratio
            set_filt1 = dict()
            for k, v in setscore.items():
                if v[0] == min_rat:
                    set_filt1[k] = v[1]

            # find sets with strongest pairs
            if len(set_filt1) > 1:
                strength = []
                for i in list(set_filt1.values()):
                    strength.append(i)
                max_str = max(strength)

                # filter by strongest pairs
                set_filt2 = []
                for k, v in set_filt1.items():
                    if v == max_str:
                        set_filt2.append(k)
                if len(set_filt2) > 1:
                    best_set = set_filt2[0]
                else:
                    best_set = set_filt2[0]

            else:
                best_set = list(set_filt1.keys())[0]

        else:
            best_set = largest_set[0]

        # dict of orthogonal pairs with data for plotting
        ort_data = dict()
        count = 0
        for p in best_set:
            ort_data[p] = pairs.get(p)
            ort_data[p].append(True)
            count = count + 1
        print(str(count) + " pairs added to dataset.")

        # find all combinations in non_pairs for orthogonality plot
        ortpep = []
        count2 = 0
        for pair in best_set:
            for pept in pair:
                ortpep.append(pept)

        for p1 in ortpep:
            for p2 in ortpep:
                if tuple([p1, p2]) in non_pairs.keys():
                    ort_data[tuple([p1, p2])] = non_pairs.get(tuple([p1, p2]))
                    ort_data[tuple([p1, p2])].append(False)
                    count2 = count2 + 1
        print(str(count2) + " peptide combinations added to dataset.")

    #######################
    ### DATA FOR EXPORT ###
    #######################

        # table with best set
        n = 0
        rows = []
        for k, v in ort_data.items():
            n = n + 1
            r = {"Peptide1": k[0],
                 "Peptide2": k[1],
                 "Kd": v[0],
                 "Interaction_score": v[1],
                 "Pair": v[2]
                 }
            rows.append(r)

        # write to disk
        table = pd.DataFrame(rows, index=None, columns=[
                                                        "Peptide1",
                                                        "Peptide2",
                                                        "Kd",
                                                        "Interaction_score",
                                                        "Pair"])

        #
        set_dir = os.path.join(dir_all, "Sets", str(orth_cutoff))
        os.makedirs(set_dir, exist_ok=True)
        fname = os.path.join(set_dir, "CC_%d" % num + str(cc_set) + "_"
                             + str(orth_cutoff) + "_orth-set-v%s.csv" % model)
        table.to_csv(fname, header=True, index=None)
        print(str(len(table)) + " total combinations written.")
        print("Written %d rows to %s" % (len(table), fname))

        # export peptide sequences
        n1 = 0
        rowsseq = []
        for pe in ortpep:
            n1 = n1 + 1
            r1 = {"Peptide": pe, "Sequence": seqs[pe]}
            rowsseq.append(r1)

        # write
        fname = os.path.join(set_dir, "CC_%d" % num + str(cc_set) + "_"
                             + str(orth_cutoff) + "_best-set_seqs-v%s.csv" % model)
        tp = pd.DataFrame(rowsseq, index=None, columns=["Peptide", "Sequence"])
        tp.to_csv(fname, header=True, index=None)
        print("Written %d rows to %s" % (len(tp), fname))

    # export table with screen data
    t = pd.DataFrame(rows_screen, index=None, columns=["Cutoff",
                                                       "Off-target",
                                                       "Homodimers",
                                                       "Weak_pairs",
                                                       "No_sets",
                                                       "Largest_set",
                                                       "Median_HP"])
    t.to_csv(out_screen, header=True, index=None)
    print("Written %d rows to %s" % (len(t), out_screen))
    print("--- %s seconds ---" % (time.time() - start_time))


# export table with full screen data
t_all = pd.DataFrame(rows_screen_all, index=None, columns=["Set",
                                                           "Cutoff",
                                                           "Off-target",
                                                           "Homodimers",
                                                           "Weak_pairs",
                                                           "No_sets",
                                                           "Largest_set",
                                                           "Median_HP"])
t_all.to_csv(out_screen_all, header=True, index=None)
print("Written %d rows to %s" % (len(t_all), out_screen_all))

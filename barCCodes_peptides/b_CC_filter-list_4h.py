import os
import sys
import pandas as pd

origin = sys.argv[1]
dirname = os.path.join(origin, "Helicities")

fname = os.path.join(dirname, "CC_list_all_4heptads_AQS_hel.csv")
df = pd.read_csv(fname, header=0, index_col=0)

# only keeps peptides that contain at least one heptad with N at a position AND at least one heptad with I at a position
keep = [("N" in s) and ("I" in s) for s in df["Sequence"]]
df_filt = df[keep]

# export filtered dataset
fname = os.path.join(dirname, "CC_list_all_4heptads_AQS_hel_filtered.csv")
df_filt.to_csv(fname, header=True)
print("Written %d rows to %s" % (len(df_filt), fname))

# split and export individual sets
sets = df_filt["Variant"].unique().tolist()
for st in sets:
    df_set = df_filt.loc[df_filt["Variant"] == str(st)]
    fname = os.path.join(dirname, "CC_list_4heptads_" + str(st) + "_hel_filtered.csv")
    df_set.to_csv(fname, header=True)
    print("Written %d rows to %s" % (len(df_set), fname))

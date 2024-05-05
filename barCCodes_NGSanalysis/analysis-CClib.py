##############
### Import ###
##############

import sys
import os
import pandas as pd
from fuzzywuzzy import fuzz
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import re


#############
### Files ###
#############


# import file with merged reads
sample_name = "libCC6"
fastq_file = sys.argv[1]

# output directory
analysis_out = sys.argv[2]
os.makedirs(analysis_out, exist_ok=True)


################
### Peptides ###
################

# define peptides used for assembly
peps = sys.argv[3]
p1_pep, p2_pep, p3_pep, null_pep = peps.split(",")


# dictionary with CC sequences
ccs = {"P0": "AIQQLAEAIAQLAQANAALAEANQALAY",
       "P1": "EIQALEEENAQLEQENAALEEEIAQLEY",
       "P3": "EIQQLEEEIAQLEQKNAALKEKNQALKY",
       "P5": "ENAALEEKIAQLKQKNAALKEEIQALEY",
       "P7": "EIQALEEKNAQLKQEIAALEEKNQALKY",
       "P9": "ENQALEQKNAQLKQEIAALEQEIAQLEY"
       }

# dictionary with expected library members
library_members = {
           "m0": [null_pep, null_pep, null_pep],
           "m1": [p1_pep, null_pep, null_pep],
           "m2": [null_pep, p2_pep, null_pep],
           "m3": [null_pep, null_pep, p3_pep],
           "m4": [p1_pep, p2_pep, null_pep],
           "m5": [p1_pep, null_pep, p3_pep],
           "m6": [null_pep, p2_pep, p3_pep],
           "m7": [p1_pep, p2_pep, p3_pep]
                   }

###################
### Definitions ###
###################


def extract_reads(fq):
    """
    Analyze Amplicon-EZ NGS reads for three-CC combinatorial libraries.

    :param fq: fastq file with merged reads
    :return: number of total reads, number of amplicon reads, number of reads matching length of correct assembly,
    dictionary with translated CC sequences for correct assemblies (translated_reads),
    dictionary with read length distribution for amplicon reads (lengths)
    """

    handle = open(fq, "rU")

    total_reads = 0  # total number of reads
    amplicon_reads = 0  # number of amplicon reads
    correct_length = 0  # reads matching length for correct assembly
    translated_reads = dict()
    lengths = dict()

    # parse file and find terminal nucleotide sequences
    for (title, sequence, quality) in FastqGeneralIterator(handle):
        total_reads = total_reads + 1  # total reads counter
        rc = str(Seq(sequence).reverse_complement())
        f = re.search("ACCGCCGCCACCATG", sequence)
        r = re.search("TACCCCTACGACGTG", sequence)
        f_rc = re.search("ACCGCCGCCACCATG", rc)
        r_rc = re.search("TACCCCTACGACGTG", rc)

    # 5` - 3` direction
        if f and r:
            amplicon_reads = amplicon_reads + 1  # amplicon reads counter
            f1, f2 = f.span()
            r1, r2 = r.span()
            length = r1 - f2

            # read length counter
            if length not in lengths:
                lengths[length] = 1
            else:
                lengths[length] = lengths[length] + 1

            if length == 450:
                correct_length = correct_length + 1  # correct assembly counter

                # translate ORF and add translated sequences for each read to dictionary
                seq = re.search("ACCGCCGCCACCATG(.*)TACCCCTACGACGTG", sequence)
                aa = Seq(seq.group(1)).translate().tostring()
                read = str(correct_length)
                translated_reads[read] = [aa[25:53], aa[68:96], aa[111:139]]

        # 3` - 5` direction
        if f_rc and r_rc:
            amplicon_reads = amplicon_reads + 1  # amplicon reads counter
            f3, f4 = f_rc.span()
            r3, r4 = r_rc.span()
            l_rc = r3 - f4

            # read length counter
            if l_rc not in lengths:
                lengths[l_rc] = 1
            else:
                lengths[l_rc] = lengths[l_rc] + 1

            if l_rc == 450:
                correct_length = correct_length + 1  # correct assembly counter

                # translate ORF and add translated sequences for each read to dictionary
                seq = re.search("ACCGCCGCCACCATG(.*)TACCCCTACGACGTG", rc)
                aa = Seq(seq.group(1)).translate().tostring()
                read = str(correct_length)
                translated_reads[read] = [aa[25:53], aa[68:96], aa[111:139]]

    sorted_lengths = dict(sorted(lengths.items(), key=lambda x: x[0]))

    return total_reads, amplicon_reads, correct_length, translated_reads, sorted_lengths


#####################
### Analyze reads ###
#####################


total_reads, amplicon_reads, match_length_reads, read_sequences, length_distribution = extract_reads(fastq_file)

percent_amplicon_reads = amplicon_reads * 100 / total_reads
percent_match = match_length_reads * 100 / amplicon_reads

print("Sample contains " + str(total_reads) + " reads.")
print(str(amplicon_reads) + " reads represent amplicon reads ("
      + str("%.1f" % percent_amplicon_reads) + "% of total reads).")
print()
print("The length of " + str(match_length_reads) + " reads matches the correct assembly length ("
      + str("%.1f" % percent_match) + "% of amplicon reads). ")


################################
### Read length distribution ###
################################


# Export length distribution as .csv
rows_length_distribution = list()
for k, v in length_distribution.items():
    if k > 0:
        rows_length_distribution.append({"Length": k, "Count": v})

lengths_df = pd.DataFrame(rows_length_distribution)
fname_length_distribution = os.path.join(analysis_out, "read_lengths.csv")
lengths_df.to_csv(fname_length_distribution, index=False)

print("Wrote file " + str(fname_length_distribution) + ".")

###################################
### Library member distribution ###
###################################


# Account for amplification and sequencing errors with fuzzy matching

fuzzy = 80  # translated CC sequence is at least 80% match to expected sequence
reads_to_members = dict()  # dictionary with reads matched to expected library members

for read, cc_combination in read_sequences.items():
    reads_to_members[read] = []

    # Position 1
    p1_fuzz_0 = fuzz.ratio(ccs[null_pep], cc_combination[0])
    p1_fuzz_1 = fuzz.ratio(ccs[p1_pep], cc_combination[0])
    if p1_fuzz_0 >= fuzzy:
        reads_to_members[read].append(null_pep)
    elif p1_fuzz_1 >= fuzzy:
        reads_to_members[read].append(p1_pep)
    else:
        reads_to_members[read].append("?")

    # Position 2
    p2_fuzz_0 = fuzz.ratio(ccs[null_pep], cc_combination[1])
    p2_fuzz_1 = fuzz.ratio(ccs[p2_pep], cc_combination[1])
    if p2_fuzz_0 >= fuzzy:
        reads_to_members[read].append(null_pep)
    elif p2_fuzz_1 >= fuzzy:
        reads_to_members[read].append(p2_pep)
    else:
        reads_to_members[read].append("?")

    # Position 3
    p3_fuzz_0 = fuzz.ratio(ccs[null_pep], cc_combination[2])
    p3_fuzz_1 = fuzz.ratio(ccs[p3_pep], cc_combination[2])
    if p3_fuzz_0 >= fuzzy:
        reads_to_members[read].append(null_pep)
    elif p3_fuzz_1 >= fuzzy:
        reads_to_members[read].append(p3_pep)
    else:
        reads_to_members[read].append("?")


# Match reads to expected library members
library_member_counts = {m: 0 for m in list(library_members.keys())}  # dictionary with counts for expected members

for read2, library_member in reads_to_members.items():
    if "?" not in library_member:
        for k, v in library_members.items():
            if v == library_member:
                library_member_counts[k] = library_member_counts[k] + 1

matched_reads = sum(library_member_counts.values())
percent_matched_reads = matched_reads * 100 / match_length_reads

print()
print(str(matched_reads) + " (" + str("%.1f" % percent_matched_reads) +
      "%) correctly assembled reads were matched to expected library members.")
print("Library member counts:")
print(pd.DataFrame.from_dict(library_member_counts, orient='index', columns=["Count"]))


# Export library member counts as .csv
rows_library_member_counts = list()
for member, count in library_member_counts.items():
    rows_library_member_counts.append({"Member_name": member,
                                       "Library_member": library_members[member],
                                       "Count": count})

counts_df = pd.DataFrame(rows_library_member_counts)
fname_member_counts = os.path.join(analysis_out, "library_member_counts.csv")
counts_df.to_csv(fname_member_counts, index=False)

print("Wrote file " + str(fname_member_counts) + ".")

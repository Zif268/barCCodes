
# BarCCodes: Coiled-coil peptide-based combinatorial labels for optical investigation of cellular architecture

Code to reproduce the study:

Lebar et al., Coiled-coil peptide-based combinatorial labels for optical investigation of cellular architecture.

## Installation

Install requirements into a new python 3 virtual environment in a directory of your choice (e.g. `~/py3cc`).

    python3 -m venv ~/py3cc
    source ~/py3cc/bin/activate
    pip install -r requirements.txt

## Analysis output

All analyses are rooted to a specified directory (e.g. `./output/`). 
    
    OUTPUT=./output/
    
## BarCCodes_peptides.

Code to reproduce coiled coil design and extraction of orthogonal sets:


### 1. Generate peptides and predict helical propensity using AGADIR.

Generate peptides with 4 and 5 heptads (27 sets per superset).

    python ./barCCodes_peptides/a_CC_generator_4heptads.py $OUTPUT
    python ./barCCodes_peptides/a_CC_generator_5heptads.py $OUTPUT


The results in the output directory should look like this.

    CC_list_AQS_4AAA.csv		CC_list_AQS_4SAQ.csv		CC_list_AQS_5QAS.csv
    CC_list_AQS_4AAQ.csv		CC_list_AQS_4SAS.csv		CC_list_AQS_5QQA.csv
    CC_list_AQS_4AAS.csv		CC_list_AQS_4SQA.csv		CC_list_AQS_5QQQ.csv
    CC_list_AQS_4AQA.csv		CC_list_AQS_4SQQ.csv		CC_list_AQS_5QQS.csv
    CC_list_AQS_4AQQ.csv		CC_list_AQS_4SQS.csv		CC_list_AQS_5QSA.csv
    CC_list_AQS_4AQS.csv		CC_list_AQS_4SSA.csv		CC_list_AQS_5QSQ.csv
    CC_list_AQS_4ASA.csv		CC_list_AQS_4SSQ.csv		CC_list_AQS_5QSS.csv
    CC_list_AQS_4ASQ.csv		CC_list_AQS_4SSS.csv		CC_list_AQS_5SAA.csv
    CC_list_AQS_4ASS.csv		CC_list_AQS_5AAA.csv		CC_list_AQS_5SAQ.csv
    CC_list_AQS_4QAA.csv		CC_list_AQS_5AAQ.csv		CC_list_AQS_5SAS.csv
    CC_list_AQS_4QAQ.csv		CC_list_AQS_5AAS.csv		CC_list_AQS_5SQA.csv
    CC_list_AQS_4QAS.csv		CC_list_AQS_5AQA.csv		CC_list_AQS_5SQQ.csv
    CC_list_AQS_4QQA.csv		CC_list_AQS_5AQQ.csv		CC_list_AQS_5SQS.csv
    CC_list_AQS_4QQQ.csv		CC_list_AQS_5AQS.csv		CC_list_AQS_5SSA.csv
    CC_list_AQS_4QQS.csv		CC_list_AQS_5ASA.csv		CC_list_AQS_5SSQ.csv
    CC_list_AQS_4QSA.csv		CC_list_AQS_5ASQ.csv		CC_list_AQS_5SSS.csv
    CC_list_AQS_4QSQ.csv		CC_list_AQS_5ASS.csv		CC_list_all_4heptads_AQS.csv
    CC_list_AQS_4QSS.csv		CC_list_AQS_5QAA.csv		CC_list_all_5heptads_AQS.csv
    CC_list_AQS_4SAA.csv		CC_list_AQS_5QAQ.csv

The next step is performed using the external tool <a href="http://agadir.crg.es">AGADIR</a>, 
for peptides in files
    
    ./output/CC_list_all_4heptads_AQS.csv
    ./output/CC_list_all_5heptads_AQS.csv

AGADIR parameters:

    pH                   7.4
    Temperature          310.15
    Ionic Strength       0.15
    
    Nterm                free
    Cterm                free

The resulting helical propensities are appended as a column `Helicity`. The resulting files are stored as:
    
    ./output/Helicities/CC_list_all_4heptads_AQS_hel.csv
    ./output/Helicities/CC_list_all_5heptads_AQS_hel.csv

The files are provided with this repository for convenience.

### 2. Filter peptide list

Input files (provided with this repository):

    ./output/Helicities/CC_list_all_4heptads_AQS_hel.csv
    ./output/Helicities/CC_list_all_5heptads_AQS_hel.csv

Filter lists and export separate files for individual sets.

    OUTPUT=./output
    python ./barCCodes_peptides/b_CC_filter-list_4h.py $OUTPUT
    python ./barCCodes_peptides/b_CC_filter-list_5h.py $OUTPUT
    

### 3. Calculate interaction scores

Calculate interaction scores for filtered peptides for models 1.0 (initial version) and 2.1 (modified version).

    # 4 heptads
    python ./barCCodes_peptides/c_CC_interaction-scores_4heptads-AQS_v1.0.py $OUTPUT
    python ./barCCodes_peptides/c_CC_interaction-scores_4heptads-AQS_v2.1.py $OUTPUT

    # 5 heptads
    python ./barCCodes_peptides/c_CC_interaction-scores_5heptads-AQS_v2.1.py $OUTPUT    


### 4. Predict orthogonality

Predict orthogonality for the given heptad number and model version.

    # 4 heptads
    python ./barCCodes_peptides/d_CC_orthogonality.py $OUTPUT 4 1.0
    python ./barCCodes_peptides/d_CC_orthogonality.py $OUTPUT 4 2.1

    # 5 heptads
    python ./barCCodes_peptides/d_CC_orthogonality.py $OUTPUT 5 2.1
    

## barCCodes_NGSanalysis

Code to replicate analysis of Amplicon EZ NGS sequencing for three-peptide combinatorial libraries.

Merged read fastq files for libCC6 and libCC7 are provided in Supplementary Data.

### 1. Download the fastq files and place them into a local folder. 

Input and output directories should be renamed accordingly.

Sample libCC6:

    OUTPUT=./output/libCC6/
    INPUT=./barCCodes_supplementary-data_SuppFig21a_libCC6-merged-reads.fastq
    
Sample libCC7:

    OUTPUT=./output/libCC7/
    INPUT=./barCCodes_supplementary-data_Fig4e_libCC7-merged-reads.fastq
    
### 2. Define peptides used for library assembly

Define peptides used for library assembly as a comma-separated list. 
First three peptides represent the building blocks and the fourth peptide represents the null peptide.
Order is important and defines the assigned positions in the peptide chain for each building block.

Sample libCC6:

    PEPS="P1,P3,P9,P7"
    
Sample libCC7:

    PEPS="P1,P3,P9,P0"
    
### 3. Run analysis
    
    python ./barCCodes_NGSanalysis/analysis-CClib.py "$INPUT" "$OUTPUT" "$PEPS"

The code writes two .csv files for each sample:

    ./output/libCC6/library_member_counts.csv
    ./output/libCC6/read_lengths.csv
    
    ./output/libCC7/library_member_counts.csv
    ./output/libCC7/read_lengths.csv
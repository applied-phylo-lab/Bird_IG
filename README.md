# Bird_IG

Proposal can be found at: https://docs.google.com/document/d/1fQ5YY_o3Em4FCX1qUgj8X3SZHkFUl0uUHCLTdSHpsR0/edit?usp=sharing

## Step-wise guide
Input: Folder with Bird data, order-subfolder structure


### Step 1: create summary_features.csv & IGH_table.tsv
``create_summary_features.R``
``create_summary_tables_clean.R``
``overview_features.R`` filter summary tables by which ones can be found in VGP tree & data overview

### Step 2: filter & clean gene files 
``filter_genes.py``

### Step 3: self-align IGH locus & create bed file from clean gene files
``IGH_self_alignment_bed.py``

### Step 4: summarize inverions
``summarize_inversions.py``

### find paralogs
- all inversions: ``find_all_inversion_paralogs.py`` Recommended
- diagonal only inversions: ``find_inversion_paralogs.py``


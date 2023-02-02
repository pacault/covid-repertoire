# Code for iScience manuscript ISCIENCE-D-22-03671R1 manuscript

`changeoToClonotypes.py` : script to convert the changeO clonotype format file to our clonotype format file

Options :

- `-i` `--input` str. output of the changeo DefineClones.py command
- `-o` `--out` str. output file name


`findcovabdab.py`: find clonotypes shared between oru repertoires and the Covabdab database

- `-c` `--covabdabfile` str. Covabdab database file (TSV format)
- `-s` `--samples` str. Samples file TSV file with columns : group name id day imgt_dir.

Examples:
```
Patient_RNA AA01 patient01 J7 /path/to/patient01/imgtdir/
Control_DNA cc01 Control01 NA /path/to/control01/imgtdir/
```

`write_clonotypes.py`: search commons clonotypes between our repertoires, and write a file with all the clonotypes.

- `-d` `--max-dif` int. Minimum mismatch for 10 aa allowed between 2 clonotypes. Default : 1
- `-n` `--min-patient-number` int. Minimum patient number that share a clonotype. Default : 4
- `-s` `--samples` str. Samples file TSV file with columns : group name id day imgt_dir.

Examples:
```
Patient_RNA AA01 patient01 J7 /path/to/patient01/imgtdir/
Control_DNA cc01 Control01 NA /path/to/control01/imgt_dir/
```

`heatmap_ms.R` : R script to plot the VJ heatmap for ARDS/HV cohorts

input: clonotypes_totaux_covabdab_fix_<ADN|ARN>.tsv

method:

  - choose ADN or ARN (beware: nonARN are not the same pts as nonADN!)
  - heatmap with library(ComplexHeatmap) 
      -> row and col annotations + gap between clusters
  - DE analysis using library(edgeR)
  - heatmap
      pre-filtering using Wald DE genes betw ARDS J0+J7 vs Control 
        and heatmap on DE VJ for all data

output:

Repertoire_<ADN|ARN>_<cond1_vs_cond2>_<allData>_heatmap.pdf
Repertoire_<ADN|ARN>_<cond1_vs_cond2>_<allData>_Top.xlsx

`covab_ms.R` : extract covabad human sequences from input file and get occurences

input: Covabdab database file (TSV format)

output: table for each VJ pair its frequency (occurence) in input file
# GenoPrimer
Automated primer design for genotyping CRISPR edited cells

## Features
- Two modes: MiSeq (250bp amplicon) and PacBio (~3000bp amplicon)
- Invokes primer3 to perform thermodynamics calculation
- Automatically relax the criteria if no primers were found.

## Inputs

- A csv file containing minimumlly four columns (with the exact names) describing each gene/primer design:
  (these)
  - Ensemble_ID (The transcript ID, e.g., ENST00000263736)  
  - Ensemble_ref (The genome/build version, e.g., GRCh38)  
  - Ensemble_chr (e.g. 2)  
  - gRNACut_in_chr (Center position of the amplicon, in the form of coordinates on the chromosome, e.g. 45389323)   
    
## Outputs:
- A csv file with the input information + new columns: 
  -  Up to three pairs of primers for each gene/row, including Tm and expected product size.
  -  A numeric number indicating how many rounds of criteria relaxation before yielding primers (column "Rounds_relax_of_primer_criteria")


&nbsp;
## Usage:
```
python GenoPrimer.py --csv test_data/test.csv
```


## Dependencies
Python 3.5-3.8 (as required by primer3-py)  
See "requirements.txt" for a complete list of dependencies



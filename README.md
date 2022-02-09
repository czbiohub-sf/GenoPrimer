# GenoPrimer
Automated primer design for genotyping CRISPR edited cells

## Features
- Two modes: MiSeq (250bp amplicon) and PacBio (~3000bp amplicon)
- Invokes primer3 to perform thermodynamics calculation
- Automatically relax the criteria if no primers were found.

## Inputs

- A csv file containing minimumlly four columns (with the exact names) describing each gene/primer design:
  (these)
  - Ensemble_ID  
    The transcript ID 
  - Ensemble_ref  
    The genome/build version
  - Ensemble_chr  
    The chromosome
  - gRNACut_in_chr  
    Center position of the amplicon, in coordinates on the chromosome

## Dependencies
Python 3.5-3.8 (as required by primer3-py)  
See "requirements.txt" for a complete list of dependencies




## Outputs:
- Genotype frequencies for each sample (.csv file)
- Alleles frequencies table (A folder containing a table of read-to-genotype assignments for each sample)

&nbsp;
## Usage:
```
python DeepGenotype.py --path2csv test_dir/test.csv --path2workDir test_dir/ --path2fastqDir test_dir/fastq_dir/
```

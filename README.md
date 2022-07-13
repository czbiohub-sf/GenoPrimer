# GenoPrimer
Automated primer design for genotyping CRISPR edited cells via amplicon sequencing

## Features
- Two modes: MiSeq (250 bp amplicon) and PacBio (~3000 bp amplicon)
- Automatically relaxes the criteria if no primers are found initially
- Invokes primer3 to perform thermodynamics calculation
- Use Blast to check unintended PCR products 
  - Autodetects system and use matching Blast executable: Linux, MacOS, Windows
  - Autodetects CPU number and multi-threads Blast search ( saves 2 CPUs for the user)

## Inputs

- A csv file containing minimumlly four columns (with the exact names) describing each gene/primer design:
  - Ensemble_ID (The transcript ID, e.g., ENST00000263736)  
  - Ensemble_ref (The genome/build version, e.g., GRCh38)  
  - Ensemble_chr (e.g. 2)  
  - gRNACut_in_chr (Center position of the amplicon, in the form of coordinates on the chromosome, e.g. 45389323)   

### [Helper script]
If you only have gRNA sequences but not their cutsites coordinates in the genome or the Ensemble IDs,
there is a helper script "get_gRNA_cutsite.py" that can obtain cutsite coordinates by mapping gRNA to the genome
more details coming soon...

## Outputs:
- A csv file with the input information + new columns: 
  -  Up to three pairs of primers for each gene/row, including Tm and expected product size.
  -  A numeric number indicating how many rounds of criteria relaxation before yielding primers (column "Rounds_relax_of_primer_criteria")

## Workflow (gene/primer design) 
<img width="1027" alt="image" src="https://user-images.githubusercontent.com/4129442/154752321-14e3f6c9-0a4c-435a-8c46-99d1a0893356.png">


&nbsp;
## Usage:
close the repository
```
git clone https://github.com/czbiohub/GenoPrimer.git
```
Go the repository directory, switch he branch if running branch other than master:
```
cd GenoPrimer
git checkout <branch you'd like to run>
```
Download the genome files 
```
cd BLAST_databases
wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```
make a Blast database, please note: [system-specific selections]
```
../bin/ncbi-blast-2.12.0+-x64-[linux/macosx/win64]/makeblastdb[.exe] -dbtype nucl -in Homo_sapiens.GRCh38.dna.primary_assembly.fa -out Homo_sapiens.GRCh38.dna.primary_assembly.fa.[linux/macosx/win64]
```
Create conda environment
```
conda env create -f environment.yml
```

You are ready to run GenoPrimer
```
cond activate GenoPrimer
python GenoPrimer.py --csv test_data/test.csv --type "MiSeq"
```






import platform
import os
import sys
import shutil
import urllib.request
import gzip
from subprocess import Popen
from Bio import SeqIO

def check_gzip_integrity(filepath):
    import gzip
    import os
    if not os.path.isfile(filepath):
        return False #file doesn't exist
    chunksize = 1024 * 1024
    with gzip.open(filepath) as g:
        try:
            while g.read(chunksize):
                pass
            return True
        except:
            return False #file corrupted

def check_blastDB(ref):
    #set BLAST binary dir
    BLAST_bin = ""
    exe_suffix = ""
    thisSystem = platform.system()
    if  thisSystem == "Windows":
        BLAST_bin = "bin/ncbi-blast-2.12.0+-x64-win64/bin/"
        exe_suffix = ".exe"
    elif thisSystem == "Darwin":
        BLAST_bin = "bin/ncbi-blast-2.12.0+-x64-macosx/bin/"
    elif thisSystem == "Linux":
        BLAST_bin = "bin/ncbi-blast-2.12.0+-x64-linux/bin/"
    else:
        sys.exit(f"Unknown operating system {thisSystem}")

    ###############################
    #check BLAST database (genome)#
    ###############################
    BLASTDB_DIR = "BLAST_databases"
    if not os.path.isdir(BLASTDB_DIR):
        os.mkdir(BLASTDB_DIR)

    if ref == "ensembl_GRCh38_latest":
        fa = "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
        prefix = "http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/"
    elif ref == "ensembl_GRCm39_latest":
        fa = "Mus_musculus.GRCm39.dna.primary_assembly.fa"
        prefix = "https://ftp.ensembl.org/pub/release-108/fasta/mus_musculus/dna/"
    elif ref == "ensembl_GRCz11_latest":
        fa = "Danio_rerio.GRCz11.dna.primary_assembly.fa"
        prefix = "https://ftp.ensembl.org/pub/release-108/fasta/danio_rerio/dna/"
    elif ref == "NCBI_refseq_GRCh38.p14":
        fa = "GCF_000001405.40_GRCh38.p14_genomic.fna"
        prefix = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/reference/GCF_000001405.40_GRCh38.p14/"
    else:
        sys.exit(f"invalid genome/ref:{ref}, possible values are: ensembl_GRCh38_latest, ensembl_GRCm39_latest, ensembl_GRCz11_latest and NCBI_refseq_GRCh38.p14")


    BLAST_db = fa + f".{thisSystem}" #make a separate db for each os system
    BLAST_db_path = os.path.join(BLASTDB_DIR, BLAST_db)

    fa_gz = fa + ".gz"
    fa_gz_path = os.path.join(BLASTDB_DIR, ref, fa_gz)
    fa = fa_gz.rstrip(".gz")

    url = prefix + fa_gz

    #check if Blast db exists
    if not (os.path.isfile(f"{BLAST_db_path}.nsq") and os.path.isfile(f"{BLAST_db_path}.ndb") and
           os.path.isfile(f"{BLAST_db_path}.nhr") and os.path.isfile(f"{BLAST_db_path}.nin") and
           os.path.isfile(f"{BLAST_db_path}.not") and os.path.isfile(f"{BLAST_db_path}.ntf") and
           os.path.isfile(f"{BLAST_db_path}.nto")):

        print(f"BLAST database not found, automatically generating BLAST databases (this is a one-time process )", flush=True)

        # check if the gz file exists, and also check integrity
        # download gzip file, download the file from `url` and save it locally under `file_name`:
        if not check_gzip_integrity(fa_gz_path): #check gzip file integrity
            print(f"...(re)downloading {fa_gz_path} sequence file, (this may take several minutes, if your internet connection is spotty, the download may fail and display an error mentioning an unexpected end of the stream)", flush=True)
            if os.path.isfile(fa_gz_path): #remove gzip file, b/c it could be a partial file
                os.remove(fa_gz_path)
            with urllib.request.urlopen(url) as response, open(fa_gz_path, 'wb') as out_file: #download gzip file
                shutil.copyfileobj(response, out_file)

        #unzip file
        if os.path.isfile(BLAST_db_path): #remove unzipped file, b/c it could be a partial file
            os.remove(BLAST_db_path)
        print(f"...Unzipping {fa_gz_path}", flush=True)
        with gzip.open(fa_gz_path) as f_in:
            with open(BLAST_db_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

        #make blastdb
        cmd = [f"{BLAST_bin}makeblastdb{exe_suffix}", "-dbtype", "nucl", "-in", f"{BLAST_db_path}"]
        print(f"...Making BLAST database, this is a one-time process", flush=True)
        p = Popen(cmd, cwd = os.getcwd(), universal_newlines=True)
        p.communicate()  # now wait plus that you can send commands to process

        os.remove(fa_gz_path) # save the gz file for making BLAST database
        os.remove(BLAST_db_path)
        print(f"Finished generating BLAST databases, this is a one-time process)", flush=True)
    return([BLAST_bin, exe_suffix, BLAST_db_path])
import platform
import os
import sys
import shutil
import urllib.request
import gzip
from subprocess import Popen

def check_blastDB_human():
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

    BLAST_db_fa = "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    BLAST_db = BLAST_db_fa + f".{thisSystem}" #make a separate db for each os system
    BLAST_db_path = os.path.join(BLASTDB_DIR, BLAST_db)
    fa_gz = "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    fa_gz_path = os.path.join(BLASTDB_DIR, fa_gz)
    fa = fa_gz.rstrip(".gz")

    url = "http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/" + fa_gz

    #check if Blast db exists
    if not (os.path.isfile(f"{BLAST_db_path}.nsq") and os.path.isfile(f"{BLAST_db_path}.ndb") and
           os.path.isfile(f"{BLAST_db_path}.nhr") and os.path.isfile(f"{BLAST_db_path}.nin") and
           os.path.isfile(f"{BLAST_db_path}.not") and os.path.isfile(f"{BLAST_db_path}.ntf") and
           os.path.isfile(f"{BLAST_db_path}.nto")):

        print(f"BLAST database not found\nDownloading {url}\nThis is a one-time process (this may take several minutes, if your internet connection is spotty, the download may fail and display an error mentioning an unexpected end of the stream)", flush=True)

        #download gzip file
        if os.path.isfile(fa_gz_path): #remove gzip file, b/c it could be a partial file
            os.remove(fa_gz_path)
        # Download the file from `url` and save it locally under `file_name`:
        with urllib.request.urlopen(url) as response, open(fa_gz_path, 'wb') as out_file:
            shutil.copyfileobj(response, out_file)

        #unzip file
        if os.path.isfile(BLAST_db_path): #remove unzipped file, b/c it could be a partial file
            os.remove(BLAST_db_path)
        print(f"unzipping {fa_gz_path}, this is a one-time process", flush=True)
        with gzip.open(fa_gz_path) as f_in:
            with open(BLAST_db_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

        #make blastdb
        cmd = [f"{BLAST_bin}makeblastdb{exe_suffix}", "-dbtype", "nucl", "-in", f"{BLAST_db_path}"]
        print(f"Making BLAST database, This is a one-time process", flush=True)
        p = Popen(cmd, cwd = os.getcwd(), universal_newlines=True)
        p.communicate()  # now wait plus that you can send commands to process

        os.remove(fa_gz_path)
        os.remove(BLAST_db_path)
    return([BLAST_bin, exe_suffix, BLAST_db_path])
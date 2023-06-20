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
        
def check_BowtieIdx(ref):
    #set Bowtie binary dir
    Bowtie_bin = ""
    exe_suffix = ""
    thisSystem = platform.system()
    if  thisSystem == "Windows":
        sys.exit(f"Windows is not supported, please use Linux or Mac")
    elif thisSystem == "Darwin":
        Bowtie_bin = "bin/bowtie-1.3.1-macos-x86_64/"
    elif thisSystem == "Linux":
        Bowtie_bin = "bin/bowtie-1.3.1-linux-x86_64/"
    else:
        sys.exit(f"Unknown operating system {thisSystem}")

    ###############################
    #check Bowtie genome indices #
    ###############################
    BowtieIdx_DIR = "Bowtie_indices"
    if not os.path.isdir(BowtieIdx_DIR):
        os.mkdir(BowtieIdx_DIR)

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


    BowtieIdx = fa + f".{thisSystem}" #make a separate db for each os system
    BowtieIdx_path = os.path.join(BowtieIdx_DIR, ref, BowtieIdx)

    fa_gz = fa + ".gz"
    fa_gz_path = os.path.join(BowtieIdx_DIR, ref, fa_gz)
    fa = fa_gz.rstrip(".gz")

    url = prefix + fa_gz


    #check if Blast db exists
    if not (os.path.isfile(f"{BowtieIdx_path}.1.ebwt") and os.path.isfile(f"{BowtieIdx_path}.2.ebwt") and
           os.path.isfile(f"{BowtieIdx_path}.3.ebwt") and os.path.isfile(f"{BowtieIdx_path}.4.ebwt") and
           os.path.isfile(f"{BowtieIdx_path}.rev.1.ebwt") and os.path.isfile(f"{BowtieIdx_path}.rev.2.ebwt")):

        print(f"Bowtie index not found, automatically generating it (this is a one-time process )", flush=True)

        # check if the gz file exists, and also check integrity
        # download gzip file, download the file from `url` and save it locally under `file_name`:
        if not check_gzip_integrity(fa_gz_path): #check gzip file integrity
            print(f"...(re)downloading {fa_gz_path} sequence file, (this may take several minutes, if your internet connection is spotty, the download may fail and display an error mentioning an unexpected end of the stream)", flush=True)
            if os.path.isfile(fa_gz_path): #remove gzip file, b/c it could be a partial file
                os.remove(fa_gz_path)
            if not os.path.isdir(os.path.dirname(fa_gz_path)):
                os.mkdir(os.path.dirname(fa_gz_path))
            with urllib.request.urlopen(url) as response, open(fa_gz_path, 'wb') as out_file: #download gzip file
                shutil.copyfileobj(response, out_file)

        #unzip file
        if os.path.isfile(BowtieIdx_path): #remove unzipped file, b/c it could be a partial file
            os.remove(BowtieIdx_path)
        print(f"...Unzipping {fa_gz_path}", flush=True)
        with gzip.open(fa_gz_path) as f_in:
            with open(BowtieIdx_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

        #generate Bowtie index
        cmd = [f"{Bowtie_bin}bowtie-build{exe_suffix}", "-f", "--offrate", "4", f"{BowtieIdx_path}", f"{BowtieIdx_path}"] # decrease offrate to 4, to speed up alignment for larger -k  
        print(f"...Making Bowtie index, this is a one-time process", flush=True)
        p = Popen(cmd, cwd = os.getcwd(), universal_newlines=True)
        p.communicate()  # now wait plus that you can send commands to process

        os.remove(fa_gz_path) # save the gz file for making BLAST database
        os.remove(BowtieIdx_path)
        print(f"Finished generating Bowtie index, this is a one-time process)", flush=True)
    return([Bowtie_bin, exe_suffix, BowtieIdx_path])

if __name__ == "__main__":
    check_BowtieIdx("ensembl_GRCh38_latest")
    check_BowtieIdx("ensembl_GRCm39_latest")
    check_BowtieIdx("ensembl_GRCz11_latest")
    #check_BowtieIdx("NCBI_refseq_GRCh38.p14")
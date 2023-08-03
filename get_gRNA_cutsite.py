import argparse
import sys
import linecache
import os
import pandas as pd
import csv
import datetime
from utils import *
import gc
import logging
from BLAST_utils import check_blastDB
from subprocess import Popen
import shutil
import urllib.request
import gzip
import traceback
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import pickle

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def parse_args():
    parser= MyParser(description='This script maps gRNA sequence to the genome to determine cutsite coordinates')
    parser.add_argument('--csv', default="", type=str, help='path to the gRNA csv file', metavar='')
    parser.add_argument('--noPAM', action='store_true', help='Ignore the requirement of a PAM following the protospacer sequence')
    config = parser.parse_args()
    if len(sys.argv)==1: # print help message if arguments are not valid
        parser.print_help()
        sys.exit(1)
    return config

config = vars(parse_args())

logging.setLoggerClass(ColoredLogger)
#logging.basicConfig()
log = logging.getLogger("gRNA mapper")
log.propagate = False
log.setLevel(logging.INFO) #set the level of warning displayed

#####################
##      main       ##
#####################
def main():
    try:
        #read input csv file
        df = pd.read_csv(os.path.join(config['csv']))

        must_have_cols = ["ref", "gene_name", "gRNA_protospacer"]
        if not all(col in df.columns for col in must_have_cols):
            log.error(f"The csv file does not contain all the required columns: ref, gene_name, gRNA_protospacer")
            sys.exit("Please fix the error(s) above and rerun the script")

        #load ENST_info
        print("Loading ENST_info")
        ENST_info = load_ENST_info()

        with open(os.path.join(f"{config['csv'].rstrip('csv')}coord.csv"), 'w') as outcsv, open(os.path.join(f"{config['csv'].rstrip('csv')}coord.log.txt"), 'w') as outlog:
            outcsv.write(",".join(list(df.columns) + ["mapping:Ensemble_chr","mapping:gRNACut_in_chr","mapping:ID","mapping:Gene_name", "mapping:Strands"])) #header
            outcsv.write("\n")
            starttime = datetime.datetime.now()
            cutsite_count = 0

            for index, row in df.iterrows():
                csvrow = [str(item) for item in row]
                genename = row["gene_name"]
                gRNAseq = row["gRNA_protospacer"]
                ref = row["ref"]
                log.info(f"processing gene {genename}, protospacer {gRNAseq}")

                gRNA_perf_match_df = get_gRNA_perf_match_in_genome(genename=genename, gRNAseq=gRNAseq, ref=ref)
                # outlog.write("\n")
                # outlog.write(gRNA_perf_match_df.to_string(index=False))
                # outlog.write("\n")

                match_chrs =[]
                match_sites = []
                match_IDs = []
                match_names = []
                match_strands = []

                for i,r in gRNA_perf_match_df.iterrows(): #go through all perfect matches
                    sstart = int(r["sstart"])
                    send =int(r["send"])
                    if sstart < send:
                        cutsite_in_chr = str(int(r["sstart"]) + 16)
                        match_strand = "+"
                    else:
                        cutsite_in_chr = str(int(r["sstart"]) - 17)
                        match_strand = "-"
                    match_chr = str(r["sseqid"])

                    #get matching IDs and gene names using the cutsite identified
                    list_of_ID_and_names = get_ENST_ID_from_site(ref=ref, Chr = match_chr, site = cutsite_in_chr,ENST_info = ENST_info)

                    tmp_IDs = []
                    tmp_names = []
                    for i in list_of_ID_and_names:
                        tmp_IDs.append(str(i[0]))
                        tmp_names.append(str(i[1]))

                    match_chrs.append(match_chr)
                    match_sites.append(cutsite_in_chr)
                    match_IDs.append( ";".join(tmp_IDs)) #there might be multiple matches
                    match_names.append(";".join(tmp_names)) #there might be multiple matches
                    match_strands.append(match_strand)

                if gRNA_perf_match_df.shape[0] > 1: # write to log for cases  more than one perfect match in the genome
                    outlog.write(",".join(csvrow) + "," + "multiple genomic targets found" + "," + "multiple genomic targets found")
                    outlog.write("\n")
                    outlog.write(gRNA_perf_match_df.to_string(index=False))
                    outlog.write("\n")

                #write to csv
                outcsv.write(",".join(csvrow) + "," + "|".join(match_chrs) + "," + "|".join(match_sites) + "," +"|".join(match_IDs) + "," + "|".join(match_names) + "," + "|".join(match_strands))
                outcsv.write("\n")
                    
                cutsite_count+=1

                if cutsite_count%50==0 and cutsite_count!=0:
                    endtime = datetime.datetime.now()
                    elapsed_sec = endtime - starttime
                    elapsed_min = elapsed_sec.seconds / 60
                    log.info(f"elapsed {elapsed_min:.2f} min, processed {cutsite_count} site(s)")
        endtime = datetime.datetime.now()
        elapsed_sec = endtime - starttime
        elapsed_min = elapsed_sec.seconds / 60
        log.info(f"finished in {elapsed_min:.2f} min, processed {cutsite_count} site(s)")

    except Exception  as e:
        print("Unexpected error:", str(sys.exc_info()))
        traceback.print_exc()
        print("additional information:", e)
        PrintException()

##########################
## function definitions ##
##########################
def get_gRNA_perf_match_in_genome(genename, gRNAseq, ref):
    numThread2use = max([1, os.cpu_count() - 2])  # for BLAST, use all CPUs except 2
    empty_query_file_flag = 1  # this value =1 if gRNA.fa is empty
    tmp_fa = f"{genename}_{gRNAseq}.fa"

    with open(tmp_fa, "w") as wfh:
        if config['noPAM']:
            PAM_len = 0
            wfh.write(f">{genename}_{gRNAseq}_\n{gRNAseq}\n")
            empty_query_file_flag = 0
        else:
            PAM_len = 3
            wfh.write(f">{genename}_{gRNAseq}_AGG\n{gRNAseq}AGG\n")
            wfh.write(f">{genename}_{gRNAseq}_CGG\n{gRNAseq}CGG\n")
            wfh.write(f">{genename}_{gRNAseq}_TGG\n{gRNAseq}TGG\n")
            wfh.write(f">{genename}_{gRNAseq}_GGG\n{gRNAseq}GGG\n")
            empty_query_file_flag = 0


    if empty_query_file_flag==0:
        #blast
        # check blastDB (human) and also return BLAST bin directory
        BLAST_bin, exe_suffix, BLAST_db_path = check_blastDB(ref)
        # specify BLAST db and query
        query = tmp_fa
        # start BLAST
        cmd = [f"{BLAST_bin}blastn{exe_suffix}", "-task", "blastn-short", "-query", f"{query}", "-db", f"{BLAST_db_path}", "-num_threads", f"{numThread2use}", "-perc_identity", "100", "-outfmt", "6 qseqid sseqid qstart qend sstart send pident mismatch", "-out", f"{query}.out"]
        p = Popen(cmd, universal_newlines=True)
        p.communicate()  # now wait for the process to finish
        os.remove(tmp_fa)

        # parse blast out
        df = pd.read_csv(f"{query}.out", sep="\t",names=["qseqid", "sseqid", "qstart", "qend", "sstart", "send", "pident", "mismatch"])
        # subset the BLAST out table to retain only perfect matches
        df = df[(df["qstart"] == 1) & (df["qend"] == (len(gRNAseq)+PAM_len)) & (df["pident"]== 100)]

        os.remove(f"{query}.out")
        return df
    else:
        return None

def read_pickle_files(file):
    if os.path.isfile(file):
        with open(file, 'rb') as handle:
            mydict = pickle.load(handle)
        return mydict
    else:
        sys.exit(f"Cannot open file: {file}")
        
def load_ENST_info():
    if not os.path.isfile(os.path.join("BLAST_databases",'ENST_info.pickle')):
        print("...preprossed ENST_info.pickle file not found, rebuilding ENST_info.pickle from GFF3 file")
        prefix = "http://ftp.ensembl.org/pub/release-110/gff3/homo_sapiens/" 
        gff3gz = "Homo_sapiens.GRCh38.110.gff3.gz"
        gff3gz_path = os.path.join("BLAST_databases",gff3gz)

        #check if gff3.gz_path exists
        if not (os.path.isfile(gff3gz_path)):
            url = prefix + gff3gz
            print(f"...GFF3 Gene annotation file not found\nDownloading {url}\nThis is a one-time process (this may take several minutes, if your internet connection is spotty, the download may fail and display an error mentioning an unexpected end of the stream)", flush=True)

            #download gzip file
            if os.path.isfile(gff3gz_path): #remove gzip file, b/c it could be a partial file
                os.remove(gff3gz_path)
            # Download the file from `url` and save it locally under `file_name`:
            with urllib.request.urlopen(url) as response, open(gff3gz_path, 'wb') as out_file:
                shutil.copyfileobj(response, out_file)

        #start extracting information from gff3gz file
        print("...Rebuilding ENST_info.pickle from GFF3 file")
        file = gff3gz_path
        #dicts
        ENST_info = dict() #dict of seq records
        ENST_exon_dict = dict()  # a dict of dicts
        ENST_exon_dict2 = dict() # a dict of req records
        ENST_CDS_dict = dict()
        ENST_codons_dict = dict()
        loc2exonID_dict = dict() # the same location may have different ENSE IDs, and need ENST to distinguish
        loc2posType = dict() # a dict that maps location to position types (e.g. 5UTR, exon intron, junction, 3UTR)

        #search patterns
        ENST_pattern = re.compile('(ENST.+?);')
        transcript_ID_pattern = re.compile('ID=transcript:(.+?);')
        Parent_pattern = re.compile('Parent=transcript:(.+?);')
        exon_id_pattern = re.compile('exon_id=(.+?);')
        CDS_id_pattern = re.compile('ID=CDS:(.+?);')
        rank_pattern = re.compile('rank=(.+?);')
        phase_pattern = re.compile('ensembl_phase=(.+?);')
        end_phase_pattern = re.compile('ensembl_end_phase=(.+?);')
        name_pattern = re.compile('Name=(.+?);')
        biotype_pattern = re.compile('biotype=(.+?);')
        version_pattern = re.compile('version=(.+?)')

        line_count = 0

        #go through gff3 file, store all ENST IDs
        print("...processing ENST IDs")
        with gzip.open(file, "rt") as fh:
            for line in fh:
                fields = line.rstrip().split("\t")
                if len(fields)>=2:
                    chr = fields[0]
                    type = fields[2]
                    m = re.search(transcript_ID_pattern, line) #require ID=ENSTXXXX
                    if m:
                        ENST_id=m.group(1)
                        if ENST_id in ENST_info:
                            sys.exit(f"{ENST_id} is already in ENST_info_dict") # ID=transcript:ENSTXXXX should be unique
                        #get transcript parent gene and it's name
                        biotype = ""
                        version = "0"
                        name = ""
                        m2 = re.search(biotype_pattern, line)
                        m3 = re.search(version_pattern, line)
                        m4 = re.search(name_pattern,line)
                        if m2: biotype = m2.group(1)
                        if m3: version = m3.group(1)
                        if m4: name = m4.group(1)
                        description = [f"{ENST_id}.{version}",biotype]
                        ENST_info[ENST_id] = SeqRecord("", id=ENST_id, description="|".join(description), name = name)
                        line_count +=1

        #go through gff3 file, and store all exons in ENST_exon_dict
        print("...processing exons")
        with gzip.open(file, "rt") as fh:
            parent_ENST_id = ""
            for line in fh:
                # print(line.rstrip())
                fields = line.rstrip().split("\t")
                if len(fields) >= 2:
                    chr = fields[0]
                    type = fields[2]
                    m2 = re.search(Parent_pattern, line)
                    if m2 and type == "exon":
                        parent_ENST_id = m2.group(1)
                        if not parent_ENST_id in ENST_info.keys():
                            sys.exit(f"{parent_ENST_id} is not found in ENST_info.keys(), offending line: {line.rstrip()} ")
                        exon_id = exon_loc = f"{parent_ENST_id}__{chr}_{fields[3]}-{fields[4]}_{fields[6]}"  # chr_start-end_strand
                        m3 = re.search(exon_id_pattern, line)
                        if m3:  # there is an exon ID
                            exon_id = m3.group(1)
                            loc2exonID_dict[exon_loc] = exon_id
                        exon_start = int(fields[3])
                        exon_end = int(fields[4])
                        exon_strand = plus_minus_strand_to_numeric(fields[6])
                        #append exon to seq record
                        ENST_info[parent_ENST_id].features.append(SeqFeature(location=FeatureLocation(exon_start,exon_end,strand=exon_strand, ref=chr),type=type, id=exon_id))

        # parse GFF3 again, extracting CDS info, and referencing exon info
        print("...processing cds")
        with gzip.open(file, "rt") as fh:
            parent_ENST_id = ""
            for line in fh:
                # print(line.rstrip())
                fields = line.rstrip().split("\t")
                if len(fields)>=2:
                    chr = fields[0]
                    type = fields[2]
                    m2 = re.search(Parent_pattern, line)
                    if m2 and type == "CDS":
                        parent_ENST_id = m2.group(1)
                        if not parent_ENST_id in ENST_info.keys():
                            sys.exit(f"{parent_ENST_id} is not found in ENST_info.keys(), offending line: {line.rstrip()} ")

                        CDS_id = f"{parent_ENST_id}__{chr}_{fields[3]}-{fields[4]}_{fields[6]}"  # chr_start-end_strand
                        # determine if CDS superimposes with an exon
                        if CDS_id in loc2exonID_dict.keys():
                            exon_id = loc2exonID_dict[CDS_id]
                            CDS_id = exon_id

                        # populate CDS info
                        CDS_start = int(fields[3])
                        CDS_end = int(fields[4])
                        CDS_strand = plus_minus_strand_to_numeric(fields[6])
                        CDS_phase = fields[7]
                        #append exon to seq record
                        #print(f"{CDS_start} {CDS_end} {CDS_strand} {type} {CDS_id}")
                        ENST_info[parent_ENST_id].features.append(SeqFeature(location=FeatureLocation(CDS_start,CDS_end,strand=CDS_strand, ref=chr),type=type, id=CDS_id))
                        # copy over additional info from exon dict
                        # if parent_ENST_id in ENST_exon_dict.keys():
                        #     if exon_id in ENST_exon_dict[parent_ENST_id].keys():
                        #         ENST_CDS_dict[parent_ENST_id][CDS_id]["exon_rank"] = ENST_exon_dict[parent_ENST_id][exon_id]["rank"]
                        #         ENST_CDS_dict[parent_ENST_id][CDS_id]["exon_phase"] = ENST_exon_dict[parent_ENST_id][exon_id]["phase"]
                        #         ENST_CDS_dict[parent_ENST_id][CDS_id]["exon_end_phase"] = ENST_exon_dict[parent_ENST_id][exon_id]["end_phase"]

        #go through ENST_info, calculate span_start span_end for each ID
        for ID in ENST_info:
            if len(ENST_info[ID].features) > 0:
                coords = [i.location.nofuzzy_end for i in ENST_info[ID].features] + [i.location.nofuzzy_start for i in ENST_info[ID].features]
                ENST_info[ID].span_start = min(coords)
                ENST_info[ID].span_end = max(coords)
                ENST_info[ID].chr = ENST_info[ID].features[0].ref
                    
        # write dict to file
        with open(os.path.join("BLAST_databases",'ENST_info.pickle'), 'wb') as handle:
            pickle.dump(ENST_info, handle, protocol=pickle.HIGHEST_PROTOCOL)
        print("...finished building ENST_info.pickle")
    
    #load ENST_info
    print("loading ENST gene model info")
    ENST_info = read_pickle_files(os.path.join("BLAST_databases",'ENST_info.pickle'))

    return ENST_info

def get_ENST_ID_from_site(ref,Chr,site,ENST_info):
    """
    return a list of lists [[ID1,name1],[ID2,name2],...]
    """
    if ref == "ensembl_GRCh38_latest":
        match_list = []
        for ID in ENST_info:
            if ENST_info[ID].chr == Chr:
                span_start = ENST_info[ID].span_start
                span_end = ENST_info[ID].span_end
                if int(span_start) <= int(site) <= int(span_end):
                    match_list.append([ID,ENST_info[ID].name])
        return match_list
    else:
        return []

def plus_minus_strand_to_numeric(input):
    if input == "+":
        return(1)
    elif input == "-":
        return(-1)
    else:
        return(0)

def PrintException():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj))


if __name__ == "__main__": main()

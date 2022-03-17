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
from BLAST_utils import check_blastDB_human
from subprocess import Popen

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def parse_args():
    parser= MyParser(description='This script maps gRNA sequence to the genome to determine cutsite coordinates')
    parser.add_argument('--csv', default="", type=str, help='path to the gRNA csv file', metavar='')
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

        with open(os.path.join(f"{config['csv'].rstrip('csv')}coord.csv"), 'w') as outcsv, open(os.path.join(f"{config['csv'].rstrip('csv')}coord.log.txt"), 'w') as outlog:
            outcsv.write(",".join(list(df.columns) + ["Ensemble_chr","gRNACut_in_chr","Ensemble_ref"])) #header
            outcsv.write("\n")
            starttime = datetime.datetime.now()
            cutsite_count = 0

            for index, row in df.iterrows():
                csvrow = [str(item) for item in row]
                genename = row["gene_name"]
                gRNAseq = row["gRNA_protospacer"]
                log.info(f"processing gene {genename}, protospacer {gRNAseq}")

                gRNA_perf_match_df = get_gRNA_perf_match_in_genome(genename=genename, gRNAseq=gRNAseq)
                # outlog.write("\n")
                # outlog.write(gRNA_perf_match_df.to_string(index=False))
                # outlog.write("\n")

                if gRNA_perf_match_df.shape[0] > 1: # more than one perfect match in the genome
                    outcsv.write(",".join(csvrow) + "," + "multiple genomic targets found" + "," + "multiple genomic targets found")
                    outcsv.write("\n")
                    outlog.write(",".join(csvrow) + "," + "multiple genomic targets found" + "," + "multiple genomic targets found")
                    outlog.write("\n")
                    outlog.write(gRNA_perf_match_df.to_string(index=False))
                    outlog.write("\n")

                else:
                    sstart = int(gRNA_perf_match_df.iloc[0]["sstart"])
                    send = int(gRNA_perf_match_df.iloc[0]["send"])
                    if sstart < send:
                        cutsite_in_chr = int(gRNA_perf_match_df.iloc[0]["sstart"]) + 17
                    else:
                        cutsite_in_chr = int(gRNA_perf_match_df.iloc[0]["sstart"]) - 17

                    outcsv.write(",".join(csvrow) + "," + str(gRNA_perf_match_df.iloc[0]["sseqid"]) + "," + str(cutsite_in_chr))
                    outcsv.write("\n")
                cutsite_count+=1
        endtime = datetime.datetime.now()
        elapsed_sec = endtime - starttime
        elapsed_min = elapsed_sec.seconds / 60
        log.info(f"finished in {elapsed_min:.2f} min, processed {cutsite_count} site(s)")

    except Exception  as e:
        print("Unexpected error:", str(sys.exc_info()))
        print("additional information:", e)
        PrintException()

##########################
## function definitions ##
##########################
def get_gRNA_perf_match_in_genome(genename, gRNAseq):
    numThread2use = max([1, os.cpu_count() - 2])  # for BLAST, use all CPUs except 2
    empty_query_file_flag = 1  # this value =1 if gRNA.fa is empty
    tmp_fa = f"{genename}_{gRNAseq}.fa"

    with open(tmp_fa, "w") as wfh:
        PAM_len = 3
        wfh.write(f">{genename}_{gRNAseq}_AGG\n{gRNAseq}AGG\n")
        wfh.write(f">{genename}_{gRNAseq}_CGG\n{gRNAseq}CGG\n")
        wfh.write(f">{genename}_{gRNAseq}_TGG\n{gRNAseq}TGG\n")
        wfh.write(f">{genename}_{gRNAseq}_GGG\n{gRNAseq}GGG\n")
        empty_query_file_flag = 0

    if empty_query_file_flag==0:
        #blast
        # check blastDB (human) and also return BLAST bin directory
        BLAST_bin, exe_suffix, BLAST_db_path = check_blastDB_human()
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

def PrintException():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj))


if __name__ == "__main__": main()

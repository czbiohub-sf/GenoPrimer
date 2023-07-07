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
import traceback
import shutil

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def parse_args():
    parser= MyParser(description='GenoPrimer')
    parser.add_argument('--csv', default="", type=str, help='path to a csv file containing three columns: ref, chr and coordinate. For details see https://github.com/czbiohub-sf/GenoPrimer#inputs, for examples, see input/example.csv', metavar='')
    parser.add_argument('--type', default="short", type=str, help='amplicon size, short:300-350bp, long: 3.5kb, default is short', metavar='')
    parser.add_argument('--thread', default="4", type=str, help='auto or an integer, auto = use max-2', metavar='')
    parser.add_argument('--outdir', default="out", type=str, help='name of the output directory relative to GenoPrimer.py', metavar='')
    #parser.add_argument('--genome', default="ensembl_GRCh38_latest", type=str, help='other accepted values are: NCBI_refseq_GRCh38.p14', metavar='')
    parser.add_argument('--db', default="Ensembl", type=str, help='name of the output directory', metavar='')
    parser.add_argument('--min_dist2edit', default = 101, type=int, help='minimum distance to the edit site', metavar='')
    parser.add_argument('--oneliner_input', default="", type=str, help='ref,chr,coordinate.  Example: ensembl_GRCh38_latest,20,17482068', metavar='')
    parser.add_argument('--aligner', default="Bowtie", type=str, help='program to align primers to the genome to check non-specific amplifications, default is Bowtie. other options: BLAST', metavar='')
    config = parser.parse_args()
    if len(sys.argv)==1: # print help message if arguments are not valid
        parser.print_help()
        sys.exit(1)
    return config

#change working directory to the script directory
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

#configs
config = vars(parse_args())
outdir = config['outdir']
oneliner = config["oneliner_input"]
min_dist2center_from_user = config["min_dist2edit"]
num_primer_return = 3
num_primers_from_Primer3 = 397 + num_primer_return

#default settings for MiSeq
prod_size_lower = 250
prod_size_upper = 350
step_size = 40
min_dist2center = 100

if config['type'] == "long":
    prod_size_lower = 3300
    prod_size_upper = 3700
    step_size = 150
    min_dist2center = 1000

# set minimum distance to edit site if a non-default value is detected
if min_dist2center_from_user != 101: 
    min_dist2center = min_dist2center_from_user

# outdir
if not os.path.isdir(outdir):
    os.makedirs(outdir)

logging.setLoggerClass(ColoredLogger)
#logging.basicConfig()
log = logging.getLogger("GenoPrimer")
log.propagate = False
log.setLevel(logging.INFO) #set the level of warning displayed

#####################
##      main       ##
#####################
def main():
    try:
        if oneliner == "":
            #check csv input
            if config["csv"] is None or config["csv"] == "":
                log.error(f"need to specify an input csv file")
                sys.exit("Please fix the error(s) above and rerun the script")

            #read input csv file
            df = pd.read_csv(os.path.join(config['csv']))
            #convert genome string to be recognizable
            df['ref'] = df['ref'].apply(get_genome_string)

            must_have_cols1 = ["ref", "chr", "coordinate"]
            flag1 = all(col in df.columns for col in must_have_cols1)
            must_have_cols2 = ["ref","mapping:Ensemble_chr","mapping:gRNACut_in_chr"]
            flag2 = all(col in df.columns for col in must_have_cols2)

            if (flag1 or flag2) == False:
                log.error(f"The csv file does not contain all the required columns: [ref, chr, coordinate] or [ref, mapping:Ensemble_chr, mapping:gRNACut_in_chr]")
                sys.exit("Please fix the error(s) above and rerun the script")
            #make a copy of the input to the output folder
            if not os.path.isfile(os.path.join(outdir,"input.csv")):
                shutil.copyfile(os.path.join(config['csv']), os.path.join(outdir,"input.csv"))
        else:
            #one-liner input
            df = pd.DataFrame([[oneliner.split(',')[0], oneliner.split(',')[1], oneliner.split(',')[2]]], columns=["ref", "chr", "coordinate"])

            must_have_cols1 = ["ref", "chr", "coordinate"]
            flag1 = all(col in df.columns for col in must_have_cols1)
            must_have_cols2 = ["ref","mapping:Ensemble_chr","mapping:gRNACut_in_chr"]
            flag2 = all(col in df.columns for col in must_have_cols2)

        #make output dir
        mkdir(outdir)

        #make output csv and begin looping over input
        with open(os.path.join(outdir, f"out.csv"), 'w') as outcsv:
            outcsv.write(",".join(list(df.columns) + ["Constraints_relaxation_iterations", "Primer Pair 1 For", "Primer Pair 1 Rev", "Primer Pair 1 For tm", "Primer Pair 1 Rev tm", "Primer Pair 1 Prod Size", "Primer Pair 2 For", "Primer Pair 2 Rev", "Primer Pair 2 For tm", "Primer Pair 2 Rev tm", "Primer Pair 2 Prod Size", "Primer Pair 3 For", "Primer Pair 3 Rev", "Primer Pair 3 For tm", "Primer Pair 3 Rev tm", "Primer Pair 3 Prod Size"])) #header
            outcsv.write("\n")
            starttime = datetime.datetime.now()
            cutsite_count = 0
            primer_count = 0
            good_primer_count = 0
            cutsite_count_noprimer = 0

            with open(os.path.join(outdir,f"log.txt"), "w") as fhlog:
                #go over each cutsite
                for index, row in df.iterrows():
                    starttime1 = datetime.datetime.now()
                    if flag1:
                        ref = row["ref"]
                        Chr = row["chr"]
                        coordinate = int(row["coordinate"])
                    if flag2:
                        ref = row["ref"]
                        Chr = str(row["mapping:Ensemble_chr"])
                        coordinate = ""
                        gene_name = row["gene_name"]
                        ENST = row["ENST"]

                        #select the correct site when multi-mapping
                        if "|" in Chr:
                            selected_idx = 0
                            for idx,item in enumerate(row["mapping:Gene_name"].split("|")): # gene name based matching, not always works
                                if gene_name !="" and not pd.isna(gene_name) and f"{gene_name}-" in item:
                                    #print(f"{gene_name}- in {item}\n")
                                    selected_idx = idx
                            for idx, item in enumerate(row["mapping:ID"].split("|")): # ENST based matching, will overwite the matching based on gene name
                                if ENST != "" and not pd.isna(ENST) and ENST in item:
                                    selected_idx = idx
                            coordinate = int(row["mapping:gRNACut_in_chr"].split("|")[selected_idx])
                            Chr = Chr.split("|")[selected_idx]
                            #print(f"{Chr} {coordinate}")
                        else:
                            coordinate = int(row["mapping:gRNACut_in_chr"])

                    log.info(f"({index+1}/{len(df.index)}) Processing cutsite:  Genome:{ref}, Chr:{Chr}, cut_coordinate: {coordinate}")
                    fhlog.write(f"({index+1}/{len(df.index)}) Processing cutsite: Genome:{ref}, Chr:{Chr}, cut_coordinate: {coordinate}\n")

                    #search for precomputed primers, end current iteration if found
                    precomputed_res = search_precomputed_results(res_dir_base = "precomputed_primers", PrimerMode = config['type'], Genome = ref, Chr = Chr, Coordinate = coordinate)
                    if precomputed_res is not None:
                        #remove the first three colums and replace with the input dataframe
                        csvrow = [str(item) for item in row]
                        res2write = ",".join(csvrow) + "," + ",".join(precomputed_res.split(",")[3:])
                        outcsv.write(res2write)
                        cutsite_count += 1
                        primer_count += 3 #TODO fix this inaccurate number
                        good_primer_count += 3 #TODO fix this inaccurate number
                        log.info(f"found precomputed primers, skip calculation for this site")
                        fhlog.write(f"found precomputed primers, skip calculation for this site")
                        continue

                    #proceed to compute primers
                    #get sequence from chromosome, get (stepsize) bp extra on each side, will progressively include in considered zone if no primers were found
                    amp_st = str(int(int(coordinate) - int(prod_size_upper)/2) - step_size*4 ) # buffer zone = step_size*4 bp
                    amp_en = str(int(int(coordinate) + int(prod_size_upper)/2) + step_size*4 ) # buffer zone = step_size*4 bp
                    chr_region = get_sequence(chromosome = str(Chr), region_left = amp_st, region_right = amp_en, genome = ref, aligner = config["aligner"]) #switched from get_ensembl_sequence() to get_sequence()

                    #design primer
                    primerlist, relaxation_count, good_primer_num = get_primers(inputSeq = str(chr_region),
                                             prod_size_lower=prod_size_lower,
                                             prod_size_upper=prod_size_upper,
                                             num_return = num_primer_return,
                                             step_size = step_size,
                                             ref = ref,
                                             chr = Chr,
                                             cut_coord = coordinate,
                                             min_dist2center = min_dist2center,
                                             num_primers_from_Primer3 = num_primers_from_Primer3,
                                             thread = config["thread"],
                                             fhlog = fhlog,
                                             outdir = outdir,
                                             aligner = config["aligner"],)

                    #process primers found
                    if primerlist is None: #no primers found
                        csvrow = [str(item) for item in row]
                        outcsv.write(",".join(csvrow) + "," + str(relaxation_count) + "," + "No qualifying primer-pairs found")
                        outcsv.write("\n")
                        cutsite_count_noprimer+=1
                    else:
                        tmp_list = flatten([[i["Lseq"], i["Rseq"], str(round(i["Ltm"],2)), str(round(i["Rtm"],2)), str(i["prodSize"])] for i in primerlist])
                        csvrow = [str(item) for item in row]
                        outcsv.write(",".join(csvrow) + "," + str(relaxation_count) + "," + ",".join(tmp_list))
                        outcsv.write("\n")
                        primer_count += len(primerlist)
                        good_primer_count += good_primer_num

                    cutsite_count += 1
                    gc.collect()

                    endtime1 = datetime.datetime.now()
                    elapsed_sec = endtime1 - starttime1
                    elapsed_min = elapsed_sec.seconds / 60
                    log.info(f"elapsed {elapsed_min:.2f} min")

                    if cutsite_count%10==0 and cutsite_count!=0:
                        endtime = datetime.datetime.now()
                        elapsed_sec = endtime - starttime
                        elapsed_min = elapsed_sec.seconds / 60
                        log.info(f"elapsed {elapsed_min:.2f} min, processed {cutsite_count} site(s), {cutsite_count_noprimer} cutcite(s) failed to yield primers")

                endtime = datetime.datetime.now()
                elapsed_sec = endtime - starttime
                elapsed_min = elapsed_sec.seconds / 60
                log.info(f"finished in {elapsed_min:.2f} min, processed {cutsite_count} site(s), found {good_primer_count} primer pair(s), outputted {primer_count} primer pair(s), {cutsite_count_noprimer} cutcite(s) failed to yield primers")
                fhlog.write(f"finished in {elapsed_min:.2f} min, processed {cutsite_count} site(s), found {good_primer_count} primer pair(s), outputted {primer_count} primer pair(s), {cutsite_count_noprimer} cutcite(s) failed to yield primers\n")
                #print(f"finished in {elapsed_min:.2f} min, processed {cutsite_count} cutsite, designed {primer_count} primers")

    except Exception  as e:
        print("Unexpected error:", str(sys.exc_info()))
        traceback.print_exc()
        print("additional information:", e)
        PrintException()

##########################
## function definitions ##
##########################
def PrintException():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj))

def mkdir(mypath):
    if not os.path.exists(mypath):
        os.makedirs(mypath)

def get_genome_string(ver):
    '''
    furnish genome string with prefix and suffix, so it will get recognized by GenoPrimer
    '''
    if not (ver.startswith("ensembl_") or ver.startswith("NCBI_refseq")):
        if config["db"] == "Ensembl":
            return f"ensembl_{ver}_latest"
        if config["db"] == "NCBI":
            return f"NCBI_refseq_{ver}.p14"
    return ver

def closest_idx(lst, K):
    return min(range(len(lst)), key = lambda i: abs(lst[i]-K))

def search_precomputed_results(res_dir_base = "precomputed_primers", PrimerMode = "short", Genome = "GRCh38", Chr = "", Coordinate = ""):
    Genome = Genome.rstrip("_latest").lstrip("ensembl_")
    result_by_chr_dir = os.path.join(res_dir_base,PrimerMode,f"ensembl_{Genome}_latest",str(Chr))
    if os.path.isdir(result_by_chr_dir):
        list_of_locsDIR = os.listdir(result_by_chr_dir)
        list_of_locs = [int(i.split("_")[0]) for i in list_of_locsDIR]
        res_dir = list_of_locsDIR[closest_idx(list_of_locs,int(Coordinate))]
        offset = abs(int(res_dir.split("_")[0]) - int(Coordinate)) 
        #print(offset)
        if offset <= 20 : # allow the target site and the precomputed site to be off by 20bp
            res_file = os.path.join(result_by_chr_dir, res_dir, "out.csv")
            if os.path.isfile(res_file):
                with open(res_file, 'r') as f: # read the precomputed results
                    header = f.readline()
                    result = f.readline()
                col_count_res = len(result.split(","))
                if col_count_res >= 9: #for valid results, change the site to the target site
                    fields = result.split(",")
                    fields[2] = str(Coordinate) # change the site in the result to the target site
                    result = ",".join(fields) 
                #no need to filter out cases where primers are not found
                return result
    return None


if __name__ == "__main__": main()

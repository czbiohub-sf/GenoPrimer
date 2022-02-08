import argparse
import sys
import linecache
import os
import pandas as pd
import csv
import datetime
from utils import *


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def parse_args():
    parser= MyParser(description='This script designs primers around the gRNA cut site')
    parser.add_argument('--csv', default="", type=str, help='path to the gRNA csv file', metavar='')
    parser.add_argument('--type', default="MiSeq", type=str, help='MiSeq:300-350bp, PacBio: 3kb', metavar='')
    config = parser.parse_args()
    if len(sys.argv)==1: # print help message if arguments are not valid
        parser.print_help()
        sys.exit(1)
    return config

config = vars(parse_args())

if config['type'] == "PacBio":
    amp_len = 3500
    prod_size_lower = 2800
    prod_size_upper = 3200
    step_size = 100

else:
    amp_len = 400
    prod_size_lower = 300
    prod_size_upper = 350
    step_size = 30

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
        #check input
        if config["csv"] is None or config["csv"] == "":
            log.error(f"need to specify an input csv file")
            sys.exit("Please fix the error(s) above and rerun the script")

        df = pd.read_csv(os.path.join(config['csv']))

        if not "Ensemble_ID" in df.columns:
            log.error(f"The csv file does not contain a column named: Ensemble_ID")
            sys.exit("Please fix the error(s) above and rerun the script")

        #read input csv file
        with open(os.path.join(f"{config['csv'].rstrip('csv')}primer.csv"), 'w') as outcsv:
            outcsv.write(",".join(list(df.columns) + ["Primer Pair 1 For", "Primer Pair 1 Rev", "Primer Pair 1 For tm", "Primer Pair 1 Rev tm", "Primer Pair 1 Prod Size", "Primer Pair 2 For", "Primer Pair 2 Rev", "Primer Pair 2 For tm", "Primer Pair 2 Rev tm", "Primer Pair 2 Prod Size", "Primer Pair 3 For", "Primer Pair 3 Rev", "Primer Pair 3 For tm", "Primer Pair 3 Rev tm", "Primer Pair 3 Prod Size"])) #header
            outcsv.write("\n")
            starttime = datetime.datetime.now()
            cutsite_count = 0
            primer_count = 0

            #go over each cutsite
            for index, row in df.iterrows():
                Ensemble_ID = row["Ensemble_ID"]
                Ensemble_spp = row["Ensemble_ref"]
                Ensemble_chr = row["Ensemble_chr"]
                gRNACut_in_chr = row["gRNACut_in_chr"]

                log.info(f"Processing cutsite: EnsembleID:{Ensemble_ID}, Genome:{Ensemble_spp}, Chr:{Ensemble_chr}, cut_coordinate: {gRNACut_in_chr}")
                #print(f"Processing cutsite: EnsembleID:{Ensemble_ID}, Genome:{Ensemble_spp}, Chr:{Ensemble_chr}, cut_coordinate: {gRNACut_in_chr}")

                #get sequence from chromosome, get 150bp extra on each side, will progressively include in considered zone if no primers were found
                amp_st = str(int(int(gRNACut_in_chr) - int(amp_len)/2) - step_size*3 ) # buffer zone = step_size*3 bp
                amp_en = str(int(int(gRNACut_in_chr) + int(amp_len)/2) + step_size*3 ) # buffer zone = step_size*3 bp
                chr_region = get_ensembl_sequence(chromosome = Ensemble_chr, region_left = amp_st, region_right = amp_en, species = "human",expand=0)

                #get left/right 10kb for off target analysis
                chr_region_left10kb = get_ensembl_sequence(chromosome = Ensemble_chr, region_left = str(int(amp_st)-10000), region_right = amp_st, species = "human",expand=0)
                chr_region_right10kb = get_ensembl_sequence(chromosome = Ensemble_chr, region_left = amp_en, region_right = str(int(amp_en)+10000), species = "human",expand=0)

                #design primer
                primerlist = get_primers(inputSeq = str(chr_region),
                                         left10kb = str(chr_region_left10kb),
                                         right10kb = str(chr_region_right10kb),
                                         prod_size_lower=prod_size_lower,
                                         prod_size_upper=prod_size_upper,
                                         num_return = 3,
                                         step_size = step_size)

                #process primers found
                if primerlist is None: #no primers found
                    csvwriter.writerow(row + ["NO primers found"])
                else:
                    tmp_list = flatten([[i["Lseq"], i["Rseq"], str(round(i["Ltm"],2)), str(round(i["Rtm"],2)), str(i["prodSize"])] for i in primerlist])
                    csvrow = [str(item) for item in row]
                    outcsv.write(",".join(csvrow) + "," + ",".join(tmp_list))
                    outcsv.write("\n")
                    primer_count += len(primerlist)

                cutsite_count += 1

            endtime = datetime.datetime.now()
            elapsed_sec = endtime - starttime
            elapsed_min = elapsed_sec.seconds / 60
            log.info(f"finished in {elapsed_min:.2f} min, processed {cutsite_count} cutsite, designed {primer_count} primers")
            #print(f"finished in {elapsed_min:.2f} min, processed {cutsite_count} cutsite, designed {primer_count} primers")
    except Exception  as e:
        print("Unexpected error:", str(sys.exc_info()))
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


if __name__ == "__main__": main()    

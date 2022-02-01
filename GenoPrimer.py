import argparse
import sys
import linecache
import primer3
import os
import pandas as pd
import requests
from Bio.Seq import Seq

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

amp_len = 400
#####################
##      main       ##
#####################    
def main():
    try:
        df = pd.read_csv(os.path.join(config['csv']))
        for index, row in df.iterrows():
            Ensemble_ID = row["Ensemble_ID"]
            Ensemble_spp = row["Ensemble_ref"]
            Ensemble_chr = row["Ensemble_chr"]
            Ensemble_chr_left_idx = row["Ensemble_chr_left_idx"]
            Ensemble_chr_right_idx = row["Ensemble_chr_right_idx"]
            gRNACut_in_chr = row["gRNACut_in_chr"]


            print(f"{Ensemble_ID} {Ensemble_spp} {Ensemble_chr} {Ensemble_chr_left_idx} {Ensemble_chr_right_idx} {gRNACut_in_chr}")


            amp_st = str(int(int(gRNACut_in_chr) - int(amp_len)/2))
            amp_en = str(int(int(gRNACut_in_chr) + int(amp_len)/2))

            chr_region = get_ensembl_sequence(chromosome = Ensemble_chr, region_left = amp_st, region_right = amp_en, species = "human",expand=0)



            get_primers(inputSeq = str(chr_region), prod_size_lower=300, prod_size_upper=350, num_return = 3)

        print("done")
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

def get_ensembl_sequence(
        chromosome,
        region_left,
        region_right,
        species,
        expand=0):
    '''
    Returns genome sequence based on chromosome range. The sequence is expanded by a flat amount on both the 5 and
    3 termini.
    '''

    base_url = "http://rest.ensembl.org"
    ext = f"/sequence/region/{species}/{chromosome}:{region_left}..{region_right}:1?expand_5prime={expand};expand_3prime={expand}"
    r = requests.get(base_url + ext, headers={"Content-Type": "text/plain"})

    if not r.ok:
        r.raise_for_status()

    # sequence = Seq(r.text, IUPACUnambiguousDNA())
    sequence = Seq(r.text)
    return sequence

def get_primers(inputSeq, prod_size_lower, prod_size_upper, num_return):
    """
    :param prod_size_lower:   product size upper bound
    :param prod_size_upper:   product size lower bound
    :param num_return: int, number of primers to return

    :return:
    a list of primers
    """
    dict_primers = primer3.bindings.designPrimers(
        {
            'SEQUENCE_ID': 'inputSeq',
            'SEQUENCE_TEMPLATE': inputSeq,
            'SEQUENCE_INCLUDED_REGION': [1, (len(inputSeq)-1)]
        },
        {
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_PICK_INTERNAL_OLIGO': 1,
            'PRIMER_INTERNAL_MAX_SELF_END': 8,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 25,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 57.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_MIN_GC': 20.0,
            'PRIMER_MAX_GC': 80.0,
            'PRIMER_MAX_POLY_X': 100,
            'PRIMER_INTERNAL_MAX_POLY_X': 100,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY': 12,
            'PRIMER_MAX_SELF_END': 8,
            'PRIMER_PAIR_MAX_COMPL_ANY': 12,
            'PRIMER_PAIR_MAX_COMPL_END': 8,
            'PRIMER_PRODUCT_SIZE_RANGE': [[prod_size_lower, prod_size_upper]],
            'PRIMER_NUM_RETURN': num_return,
        })

    # go through primers found
    primer_pair_num = dict_primers['PRIMER_PAIR_NUM_RETURNED']
    for i in range(0, primer_pair_num):
        Rank = i
        Lseq = dict_primers["PRIMER_LEFT_{}_SEQUENCE".format(i)]
        Rseq = dict_primers["PRIMER_RIGHT_{}_SEQUENCE".format(i)]
        Ltm = dict_primers["PRIMER_LEFT_{}_TM".format(i)]
        Rtm = dict_primers["PRIMER_RIGHT_{}_TM".format(i)]
        prodSize = dict_primers["PRIMER_PAIR_{}_PRODUCT_SIZE".format(i)]
        print("{},{},{},{},{},{}".format(i, Lseq, Rseq, Ltm, Rtm, prodSize))

if __name__ == "__main__": main()    

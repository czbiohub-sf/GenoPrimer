import argparse
import sys
import linecache
import primer3
import os
import pandas as pd
import requests
from Bio.Seq import Seq
import re
import csv

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

else:
    amp_len = 400
    prod_size_lower = 300
    prod_size_upper = 350

#####################
##      main       ##
#####################    
def main():
    try:
        df = pd.read_csv(os.path.join(config['csv']))
        with open(os.path.join(f"{config['csv'].rstrip('csv')}primer.csv"), 'w') as outcsv:
            outcsv.write(",".join(list(df.columns) + ["Primer Pair 1 For", "Primer Pair 1 Rev", "Primer Pair 1 For tm", "Primer Pair 1 Rev tm", "Primer Pair 1 Prod Size", "Primer Pair 2 For", "Primer Pair 2 Rev", "Primer Pair 2 For tm", "Primer Pair 2 Rev tm", "Primer Pair 2 Prod Size", "Primer Pair 3 For", "Primer Pair 3 Rev", "Primer Pair 3 For tm", "Primer Pair 3 Rev tm", "Primer Pair 3 Prod Size"])) #header
            outcsv.write("\n")
            for index, row in df.iterrows():
                Ensemble_ID = row["Ensemble_ID"]
                Ensemble_spp = row["Ensemble_ref"]
                Ensemble_chr = row["Ensemble_chr"]
                Ensemble_chr_left_idx = row["Ensemble_chr_left_idx"]
                Ensemble_chr_right_idx = row["Ensemble_chr_right_idx"]
                gRNACut_in_chr = row["gRNACut_in_chr"]

                print(f"{Ensemble_ID} {Ensemble_spp} {Ensemble_chr} {Ensemble_chr_left_idx} {Ensemble_chr_right_idx} {gRNACut_in_chr}")

                amp_st = str(int(int(gRNACut_in_chr) - int(amp_len)/2) - 150) # buffer zone = 300bp
                amp_en = str(int(int(gRNACut_in_chr) + int(amp_len)/2) + 150) # buffer zone = 300bp

                chr_region = get_ensembl_sequence(chromosome = Ensemble_chr, region_left = amp_st, region_right = amp_en, species = "human",expand=0)
                chr_region_left10kb = get_ensembl_sequence(chromosome = Ensemble_chr, region_left = str(int(amp_st)-10000), region_right = amp_st, species = "human",expand=0)
                chr_region_right10kb = get_ensembl_sequence(chromosome = Ensemble_chr, region_left = amp_en, region_right = str(int(amp_en)+10000), species = "human",expand=0)

                primerlist = get_primers(inputSeq = str(chr_region),
                                         left10kb = str(chr_region_left10kb),
                                         right10kb = str(chr_region_right10kb),
                                         prod_size_lower=prod_size_lower,
                                         prod_size_upper=prod_size_upper,
                                         num_return = 3)

                if primerlist is None: #no primers found
                    csvwriter.writerow(row + ["NO primers found"])
                else:
                    tmp_list = flatten([[i["Lseq"], i["Rseq"], str(i["Ltm"]), str(i["Rtm"]), str(i["prodSize"])] for i in primerlist])
                    csvrow = [str(item) for item in row]
                    outcsv.write(",".join(csvrow) + "," + ",".join(tmp_list))
                    outcsv.write("\n")

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

def flatten(t):
    return [item for sublist in t for item in sublist]

def get_ensembl_sequence(chromosome,region_left,region_right,species,expand=0):
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

def get_primers(inputSeq, left10kb, right10kb, prod_size_lower, prod_size_upper, num_return):
    """
    :param prod_size_lower:   product size upper bound
    :param prod_size_upper:   product size lower bound
    :param num_return: int, number of primers to return

    :return:
    a list of primers
    """
    #create dicts as inputs for primer3
    thermo_dict = get_default_thermo_dict()

    User_dict1={
            'SEQUENCE_ID': 'inputSeq',
            'SEQUENCE_TEMPLATE': inputSeq,
            'SEQUENCE_INCLUDED_REGION': [1+150, (len(inputSeq)-1)-150] # 150bp buffer zone on each side
    }
    User_dict2 = {
            'PRIMER_PRODUCT_SIZE_RANGE': [[prod_size_lower, prod_size_upper]],
            'PRIMER_NUM_RETURN': num_return + 5 # get 5 extra primer pairs
    }

    dict_primers = primer3.bindings.designPrimers(
        User_dict1,
        {**thermo_dict, **User_dict2} # these two dicts needs to be merged
)
    #check unintended products
    dict_primers = check_unintended_products(dict_primers = dict_primers,
                                             left10kb = left10kb,
                                             right10kb = right10kb,
                                             len_input = User_dict1['SEQUENCE_INCLUDED_REGION'][1] - User_dict1['SEQUENCE_INCLUDED_REGION'][0])

    #get primer number
    primer_pair_num = dict_primers['PRIMER_PAIR_NUM_RETURNED'] - len(dict_primers["UNSPECIFIC_PRIMER_PAIR_idx"])



    if primer_pair_num < 1: #no primers found
        relaxation_count = 1
        while (relaxation_count < 7 or primer_pair_num < 1):
            #alternate between two relaxations
            if relaxation_count % 2 == 0:
                thermo_dict = relax_MIN_MAX_TM(thermo_dict) #relax thermodynamics
            else:
                User_dict1, User_dict2 = relax_amp_size(User_dict1, User_dict2) #relax amplicon size
            dict_primers = primer3.bindings.designPrimers(User_dict1, {**thermo_dict, **User_dict2})
            # check unintended products
            dict_primers = check_unintended_products(dict_primers=dict_primers,
                                                     left10kb=left10kb,
                                                     right10kb=right10kb,
                                                     len_input=User_dict1['SEQUENCE_INCLUDED_REGION'][1] - User_dict1['SEQUENCE_INCLUDED_REGION'][0])
            # get primer number
            primer_pair_num = dict_primers['PRIMER_PAIR_NUM_RETURNED'] - len(dict_primers["UNSPECIFIC_PRIMER_PAIR_idx"]) # check primers found using relaxed criteria
            relaxation_count += 1

    if primer_pair_num < 1:  # still no primers found
        return
    else: #primer found:
        primer_list = []
        for i in range(0, min(primer_pair_num, num_return)):
            if not i in dict_primers["UNSPECIFIC_PRIMER_PAIR_idx"]: #avoid getting unspecific primers
                tmpDict ={}
                Rank = i+1
                tmpDict["Rank"] = Rank
                tmpDict["Lseq"] = dict_primers["PRIMER_LEFT_{}_SEQUENCE".format(i)]
                tmpDict["Rseq"] = dict_primers["PRIMER_RIGHT_{}_SEQUENCE".format(i)]
                tmpDict["Ltm"] = dict_primers["PRIMER_LEFT_{}_TM".format(i)]
                tmpDict["Rtm"] = dict_primers["PRIMER_RIGHT_{}_TM".format(i)]
                tmpDict["prodSize"] = dict_primers["PRIMER_PAIR_{}_PRODUCT_SIZE".format(i)]
                primer_list.append(tmpDict)
        return primer_list

def get_default_thermo_dict():
    return {
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
        }

def relax_MIN_MAX_TM(thermo_dict):
    """
    decrease PRIMER_MIN_TM by 1
    increase PRIMER_MAX_TM by 1
    :param thermo_dict:
    :return: thermo_dict
    """
    thermo_dict['PRIMER_MIN_TM'] -= 1
    thermo_dict['PRIMER_MAX_TM'] += 1
    return thermo_dict

def relax_amp_size(User_dict1,User_dict2):
    User_dict1["SEQUENCE_INCLUDED_REGION"] = [User_dict1["SEQUENCE_INCLUDED_REGION"][0] - 50,
                                              User_dict1["SEQUENCE_INCLUDED_REGION"][1] + 50]
    User_dict2["PRIMER_PRODUCT_SIZE_RANGE"] = [[User_dict2["PRIMER_PRODUCT_SIZE_RANGE"][0][0] + 50,
                                                 User_dict2["PRIMER_PRODUCT_SIZE_RANGE"][0][1] + 50]]

    return [User_dict1, User_dict2]

def check_unintended_products(dict_primers, left10kb, right10kb, len_input):
    """
    checks unintended_products and remove primers having unintended products
    :param dict_primers:
    :param left10kb:
    :param right10kb:
    :return: dict_primers(updated) #add a key PRIMER_PAIR_NUM_RETURNED_SPECIFIC
    """
    dict_primers["UNSPECIFIC_PRIMER_PAIR_idx"] = []
    for i in range(0, dict_primers["PRIMER_PAIR_NUM_RETURNED"]):
        Lseq = dict_primers["PRIMER_LEFT_{}_SEQUENCE".format(i)]
        Rseq = dict_primers["PRIMER_RIGHT_{}_SEQUENCE".format(i)]

        flag1 = check_unintended_primer_pairing(Lseq, Rseq, left10kb, min_span=80, max_span = len_input)
        flag2 = check_unintended_primer_pairing(Lseq, Rseq, right10kb, min_span=80, max_span = len_input)

        if flag1 or flag2: # unintended product
            dict_primers["UNSPECIFIC_PRIMER_PAIR_idx"].append(i)

    return dict_primers

def revcom(seq):
    return str(Seq(seq).reverse_complement())

def check_unintended_primer_pairing(For, Rev, seq, min_span, max_span): ##check if any matches produces amplicons
    """
    :param For:
    :param Rev:
    :param seq:
    :param min_span:
    :param max_span:
    :return: True if found unintended amplification, and False if otherwise
    """
    For_pat = re.compile(For[-12:], re.IGNORECASE)
    RevRC_pat = re.compile(revcom(Rev)[0:12], re.IGNORECASE)

    For_matches = []
    RevRC_matches = []
    for m in For_pat.finditer(seq):
        For_matches.append(m.span()[0])
    for m in RevRC_pat.finditer(seq):
        RevRC_matches.append(m.span()[1])

    ## -> go through matches
    if For_matches and RevRC_matches:
        for i in For_matches:
            for j in RevRC_matches:
                if min_span <= (int(j) - int(i)) <= max_span:
                    return True

    #check matches in the revcom strand
    For_matches = []
    RevRC_matches = []
    for m in For_pat.finditer(revcom(seq)):
        For_matches.append(m.span()[0])
    for m in RevRC_pat.finditer(revcom(seq)):
        RevRC_matches.append(m.span()[1])

    ## -> go through matches
    if For_matches and RevRC_matches:
        for i in For_matches:
            for j in RevRC_matches:
                if min_span <= (int(j) - int(i)) <= max_span:
                    return True
    return False

if __name__ == "__main__": main()    

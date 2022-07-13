import primer3
import requests
from Bio.Seq import Seq
import re
import math
import csv
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
from check_unintend_prod import check_unintended_products

def flatten(t):
    return [item for sublist in t for item in sublist]

def get_ensembl_sequence(chromosome,region_left,region_right,species,expand=0):
    '''
    Returns genome sequence based on chromosome range. The sequence is expanded by a flat amount on both the 5 and
    3 termini.
    '''
    #retry strategy
    retry_strategy = Retry(
        total=7,
        backoff_factor=2, #how long the processes will sleep = {backoff factor} * (2 ** ({number of total retries} - 1));  2 seconds = 1, 2, 4, 8, 16, 32, 64, 128, 256, 512
        status_forcelist=[429, 500, 502, 503, 504],
        method_whitelist=["HEAD", "GET", "OPTIONS"]
    )
    adapter = HTTPAdapter(max_retries=retry_strategy)
    http = requests.Session()
    http.mount("https://", adapter)
    http.mount("http://", adapter)

    base_url = "http://rest.ensembl.org"
    ext = f"/sequence/region/{species}/{chromosome}:{region_left}..{region_right}:1?expand_5prime={expand};expand_3prime={expand}"
    r = http.get(base_url + ext, headers={"Content-Type": "text/plain"})

    if not r.ok:
        r.raise_for_status()

    # sequence = Seq(r.text, IUPACUnambiguousDNA())
    sequence = Seq(r.text)
    return sequence

def get_primers(inputSeq, prod_size_lower, prod_size_upper, num_return, step_size, chr, cut_coord, min_dist2center, num_primers_from_Primer3,fhlog):
    """
    :param prod_size_lower:   product size upper bound
    :param prod_size_upper:   product size lower bound
    :param num_return: int, number of primers to return

    :return:
    a list of primers
    """
    ##################
    #setup parameters#
    ##################
    min_dist2center = min_dist2center

    #create dicts as inputs for primer3
    thermo_dict = get_default_thermo_dict()

    User_dict1={
            'SEQUENCE_ID': 'inputSeq',
            'SEQUENCE_TEMPLATE': inputSeq,
            'SEQUENCE_EXCLUDED_REGION': [[1,3*step_size],
                                         [math.floor(len(inputSeq)/2) - min_dist2center, 2*min_dist2center],
                                         [len(inputSeq)-3*step_size-1, 3*step_size]]
            #'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST ': [[1 + step_size*3,                                                                    #Forward primer region start
            #                                         math.floor(len(inputSeq)/2) - min_dist2center - 1 - step_size*3,                     #Forward primer region length
            #                                         math.floor(len(inputSeq)/2) + min_dist2center,                                       #Reverse primer region start
            #                                         (len(inputSeq)-1) - step_size*3 - math.floor(len(inputSeq)/2) -  min_dist2center]]  #Reverse primer region length
    }
    User_dict2 = {
            'PRIMER_PRODUCT_SIZE_RANGE': [[max(prod_size_lower,2*min_dist2center), prod_size_upper]],  #product lower size can't be smaller than 2*min_dist2center
            'PRIMER_NUM_RETURN': num_primers_from_Primer3 # get 7 extra primer pairs
    }

    #keep a list of primers that failed specificity check, to avoid repeated work
    nonspecific_primers = {}

    ################
    #search primers#
    ################
    dict_primers = primer3.bindings.designPrimers(
        User_dict1,
        {**thermo_dict, **User_dict2}) # these two dicts needs to be merged
    #check unintended products
    dict_primers = check_unintended_products(dict_primers = dict_primers, len_input = prod_size_upper, cut_chr = chr, cut_coord = cut_coord, nonspecific_primers=nonspecific_primers, fhlog = fhlog)
    nonspecific_primers = populate_nonspecific_primer_list(dict_primers, nonspecific_primers)

    #get primer number
    primer_pair_num = dict_primers['PRIMER_PAIR_NUM_RETURNED'] - len(dict_primers["UNSPECIFIC_PRIMER_PAIR_idx"])

    # check if all primers fail because of unintended product
    if dict_primers['PRIMER_PAIR_NUM_RETURNED'] == len(dict_primers["UNSPECIFIC_PRIMER_PAIR_idx"]) and not dict_primers["PRIMER_PAIR_NUM_RETURNED"] == 0:
        print("All primers failed because of having unintended PCR products, checking 100 candidate primers for their unintended PCR products ", flush=True)
        fhlog.write("All primers failed because of having unintended PCR products, checking 100 candidate primers for their unintended PCR products \n")
        # get more primer candiates to check
        User_dict2 = primer_num_eq_800(User_dict2)
        # design primers
        dict_primers = primer3.bindings.designPrimers(User_dict1, {**thermo_dict, **User_dict2})
        # check unintended products
        dict_primers = check_unintended_products(dict_primers=dict_primers, len_input=prod_size_upper, cut_chr=chr, cut_coord=cut_coord, nonspecific_primers=nonspecific_primers, fhlog = fhlog )
        nonspecific_primers = populate_nonspecific_primer_list(dict_primers, nonspecific_primers)

        # get primer number
        primer_pair_num = dict_primers['PRIMER_PAIR_NUM_RETURNED'] - len(
            dict_primers["UNSPECIFIC_PRIMER_PAIR_idx"])  # check primers found using relaxed criteria
        User_dict2 = reset_primer_num(User_dict2, num_return)

    #########################
    #Begin relaxing criteria#
    #########################
    relaxation_count = 0
    if primer_pair_num < 1: #no primers found, start relaxation of criteria
        while (relaxation_count < 7 and primer_pair_num < 1):
            print(f"[Warning] No primers found, relaxing criteria, iteration={relaxation_count+1}", flush=True)
            fhlog.write(f"[Warning] No primers found, relaxing criteria, iteration={relaxation_count + 1}\n")
            #alternate between two relaxations
            if relaxation_count % 2 == 0:
                thermo_dict = relax_MIN_MAX_TM(thermo_dict) #relax thermodynamics
            else:
                User_dict1, User_dict2 = relax_amp_size(User_dict1, User_dict2, step_size) #relax amplicon size

            dict_primers = primer3.bindings.designPrimers(User_dict1, {**thermo_dict, **User_dict2})
            # check unintended products
            dict_primers = check_unintended_products(dict_primers=dict_primers, len_input=prod_size_upper, cut_chr = chr, cut_coord = cut_coord, nonspecific_primers=nonspecific_primers, fhlog = fhlog )
            nonspecific_primers = populate_nonspecific_primer_list(dict_primers, nonspecific_primers)

            # get primer number
            primer_pair_num = dict_primers['PRIMER_PAIR_NUM_RETURNED'] - len(dict_primers["UNSPECIFIC_PRIMER_PAIR_idx"]) # check primers found using relaxed criteria
            relaxation_count += 1

            #check if all primers fail because of unintended product
            if dict_primers['PRIMER_PAIR_NUM_RETURNED'] == len(dict_primers["UNSPECIFIC_PRIMER_PAIR_idx"]) and not dict_primers["PRIMER_PAIR_NUM_RETURNED"]==0:
                print("All primers failed because of having unintended PCR products, checking 1000 candidate primers for their unintended PCR products ", flush=True)
                fhlog.write("All primers failed because of having unintended PCR products, checking 1000 candidate primers for their unintended PCR products ")
                #get more primer candiates to check
                User_dict2 = primer_num_eq_800(User_dict2)
                #design primers
                dict_primers = primer3.bindings.designPrimers(User_dict1, {**thermo_dict, **User_dict2})
                # check unintended products
                dict_primers = check_unintended_products(dict_primers=dict_primers, len_input=prod_size_upper, cut_chr = chr, cut_coord = cut_coord, nonspecific_primers=nonspecific_primers, fhlog = fhlog )
                nonspecific_primers = populate_nonspecific_primer_list(dict_primers, nonspecific_primers)

                # get primer number
                primer_pair_num = dict_primers['PRIMER_PAIR_NUM_RETURNED'] - len(dict_primers["UNSPECIFIC_PRIMER_PAIR_idx"]) # check primers found using relaxed criteria
                User_dict2 = reset_primer_num(User_dict2, num_return)


    if primer_pair_num<1: #no primer found after a series relaxation, do a final relaxation of distance to center
        print(f"[Warning] No primers found, relaxing criteria, iteration={relaxation_count + 1}", flush=True)
        fhlog.write(f"[Warning] No primers found, relaxing criteria, iteration={relaxation_count + 1}\n")
        length_closer_tocenter = int(min_dist2center * 0.8) # for 100bp, lengh_closer_tocenter = 20
        User_dict1, User_dict2 = relax_dist2center(User_dict1 = User_dict1, User_dict2 = User_dict2, length_closer_tocenter = length_closer_tocenter)
        dict_primers = primer3.bindings.designPrimers(User_dict1, {**thermo_dict, **User_dict2})
        # check unintended products
        dict_primers = check_unintended_products(dict_primers=dict_primers, len_input=prod_size_upper, cut_chr = chr, cut_coord = cut_coord, nonspecific_primers=nonspecific_primers, fhlog = fhlog )
        nonspecific_primers = populate_nonspecific_primer_list(dict_primers, nonspecific_primers)

        # get primer number
        primer_pair_num = dict_primers['PRIMER_PAIR_NUM_RETURNED'] - len(
            dict_primers["UNSPECIFIC_PRIMER_PAIR_idx"])  # check primers found using relaxed criteria
        relaxation_count += 1

    ########
    #output#
    ########
    if primer_pair_num < 1:  # still no primers found
        return [None, relaxation_count, 0]
    else: #primer found:
        primer_list = []
        for i in range(0, dict_primers['PRIMER_PAIR_NUM_RETURNED']): # go through all primers
            if not i in dict_primers["UNSPECIFIC_PRIMER_PAIR_idx"]: #avoid getting unspecific primers
                if len(primer_list) == num_return: # return list if full, check if current primer is further away from the 1st primer, if so, replace the last primer
                    new_L_st = dict_primers["PRIMER_LEFT_{}".format(i)][0]
                    new_R_st = dict_primers["PRIMER_RIGHT_{}".format(i)][0]
                    first_L_st = primer_list[0]["Lstart"]
                    first_R_st = primer_list[0]["Rstart"]
                    last_L_st = primer_list[len(primer_list)-1]["Lstart"]
                    last_R_st = primer_list[len(primer_list)-1]["Rstart"]

                    new_dist     = abs(new_L_st - first_L_st) + abs(new_R_st - first_R_st) # new dist between primer pairs
                    current_dist = abs(last_L_st - first_L_st) + abs(last_R_st - first_R_st) # current dist between primer pairs

                    if new_dist > current_dist: # replace the last primer pair with new one if the new dist is larger
                        tmpDict = {}
                        Rank = i + 1
                        tmpDict["Rank"] = Rank
                        tmpDict["Lseq"] = dict_primers["PRIMER_LEFT_{}_SEQUENCE".format(i)]
                        tmpDict["Rseq"] = dict_primers["PRIMER_RIGHT_{}_SEQUENCE".format(i)]
                        tmpDict["Ltm"] = dict_primers["PRIMER_LEFT_{}_TM".format(i)]
                        tmpDict["Rtm"] = dict_primers["PRIMER_RIGHT_{}_TM".format(i)]
                        tmpDict["prodSize"] = dict_primers["PRIMER_PAIR_{}_PRODUCT_SIZE".format(i)]
                        tmpDict["Lstart"] = dict_primers["PRIMER_LEFT_{}".format(i)][0]
                        tmpDict["Rstart"] = dict_primers["PRIMER_RIGHT_{}".format(i)][0]

                        primer_list[len(primer_list) - 1] = tmpDict # replace the last primer pair with current one

                else: # add primer to return list
                    tmpDict ={}
                    Rank = i+1
                    tmpDict["Rank"] = Rank
                    tmpDict["Lseq"] = dict_primers["PRIMER_LEFT_{}_SEQUENCE".format(i)]
                    tmpDict["Rseq"] = dict_primers["PRIMER_RIGHT_{}_SEQUENCE".format(i)]
                    tmpDict["Ltm"] = dict_primers["PRIMER_LEFT_{}_TM".format(i)]
                    tmpDict["Rtm"] = dict_primers["PRIMER_RIGHT_{}_TM".format(i)]
                    tmpDict["prodSize"] = dict_primers["PRIMER_PAIR_{}_PRODUCT_SIZE".format(i)]
                    tmpDict["Lstart"] = dict_primers["PRIMER_LEFT_{}".format(i)][0]
                    tmpDict["Rstart"] = dict_primers["PRIMER_RIGHT_{}".format(i)][0]
                    primer_list.append(tmpDict)
        good_primer_num = dict_primers['PRIMER_PAIR_NUM_RETURNED'] - len(dict_primers["UNSPECIFIC_PRIMER_PAIR_idx"])
        return [primer_list, relaxation_count, good_primer_num]

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
    thermo_dict['PRIMER_MIN_TM'] -= 0.5
    thermo_dict['PRIMER_MAX_TM'] += 0.5
    return thermo_dict

def relax_amp_size(User_dict1,User_dict2,step_size):
    #User_dict1["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"] = [[User_dict1["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"][0][0] - step_size,
    #                                                      User_dict1["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"][0][1] + step_size,   # start - stepsize,  len + stepsize
    #                                                      User_dict1["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"][0][2],
    #                                                      User_dict1["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"][0][3]+step_size]]  #start, len+stepsize
    User_dict1['SEQUENCE_EXCLUDED_REGION'][0][1] = User_dict1['SEQUENCE_EXCLUDED_REGION'][0][1] - step_size # left len - stepsize
    User_dict1['SEQUENCE_EXCLUDED_REGION'][2][0] = User_dict1['SEQUENCE_EXCLUDED_REGION'][2][0] + step_size # right start + stepsize
    User_dict1['SEQUENCE_EXCLUDED_REGION'][2][1] = User_dict1['SEQUENCE_EXCLUDED_REGION'][2][1] - step_size # right len - stepsize

    User_dict2["PRIMER_PRODUCT_SIZE_RANGE"] = [[User_dict2["PRIMER_PRODUCT_SIZE_RANGE"][0][0], #product size lower stay the same
                                                 User_dict2["PRIMER_PRODUCT_SIZE_RANGE"][0][1] + step_size*2]] #product size upper increases
    return [User_dict1, User_dict2]

def relax_dist2center(User_dict1,User_dict2, length_closer_tocenter):
    #User_dict1["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"] = [[User_dict1["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"][0][0],
    #                                                      User_dict1["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"][0][1] + length_closer_tocenter, # start ,  len + length_closer_tocenter
    #                                                      User_dict1["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"][0][2] - length_closer_tocenter,
    #                                                      User_dict1["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"][0][3] + length_closer_tocenter]]  #start - length_closer_tocenter, len + length_closer_tocenter
    User_dict1['SEQUENCE_EXCLUDED_REGION'][1][0] = User_dict1['SEQUENCE_EXCLUDED_REGION'][1][0] + length_closer_tocenter # middle start + length_closer_tocenter
    User_dict1['SEQUENCE_EXCLUDED_REGION'][1][1] = User_dict1['SEQUENCE_EXCLUDED_REGION'][1][1] - 2*length_closer_tocenter # middle len - 2*length_closer_tocenter

    User_dict2["PRIMER_PRODUCT_SIZE_RANGE"] = [[User_dict2["PRIMER_PRODUCT_SIZE_RANGE"][0][0] - 2*length_closer_tocenter, #product size lower decrease
                                                 User_dict2["PRIMER_PRODUCT_SIZE_RANGE"][0][1] ]] #product size stay the same

    return [User_dict1, User_dict2]

def primer_num_eq_1000(User_dict2):
    User_dict2["PRIMER_NUM_RETURN"] = 1000
    return User_dict2

def primer_num_eq_800(User_dict2):
    User_dict2["PRIMER_NUM_RETURN"] = 800
    return User_dict2

def reset_primer_num(User_dict2,num_return):
    User_dict2["PRIMER_NUM_RETURN"] = num_return + 5
    return User_dict2

def revcom(seq):
    return str(Seq(seq).reverse_complement())

def check_unintended_primer_pairing(For, Rev, min_span, max_span): ##check if any matches produces amplicons
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

def populate_nonspecific_primer_list(dict_primers, nonspecific_primers):
    for idx in range(0, dict_primers['PRIMER_PAIR_NUM_RETURNED']):
        if str(idx) in dict_primers['UNSPECIFIC_PRIMER_PAIR_idx']:
            Lseq = dict_primers["PRIMER_LEFT_{}_SEQUENCE".format(idx)]
            Rseq = dict_primers["PRIMER_RIGHT_{}_SEQUENCE".format(idx)]
            nonspecific_primers[f"{Lseq}{Rseq}"]=1
    return nonspecific_primers

#################
#custom logging #
#################
import logging

BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)
#The background is set with 40 plus the number of the color, and the foreground with 30

#These are the sequences need to get colored ouput
RESET_SEQ = "\033[0m"
COLOR_SEQ = "\033[1;%dm"
BOLD_SEQ = "\033[1m"

def formatter_message(message, use_color = True):
    if use_color:
        message = message.replace("$RESET", RESET_SEQ).replace("$BOLD", BOLD_SEQ)
    else:
        message = message.replace("$RESET", "").replace("$BOLD", "")
    return message

COLORS = {
    'WARNING': YELLOW,
    'INFO': WHITE,
    'DEBUG': BLUE,
    'CRITICAL': YELLOW,
    'ERROR': RED
}

class ColoredFormatter(logging.Formatter):
    def __init__(self, msg, use_color = True):
        logging.Formatter.__init__(self, msg)
        self.use_color = use_color

    def format(self, record):
        levelname = record.levelname
        if self.use_color and levelname in COLORS:
            levelname_color = COLOR_SEQ % (30 + COLORS[levelname]) + levelname + RESET_SEQ
            record.levelname = levelname_color
        return logging.Formatter.format(self, record)

# Custom logger class with multiple destinations
class ColoredLogger(logging.Logger):
    FORMAT = "[$BOLD%(name)-1s$RESET][%(levelname)-1s]  %(message)s " #($BOLD%(filename)s$RESET:%(lineno)d)
    COLOR_FORMAT = formatter_message(FORMAT, True)
    def __init__(self, name):
        logging.Logger.__init__(self, name, logging.DEBUG)

        color_formatter = ColoredFormatter(self.COLOR_FORMAT)

        console = logging.StreamHandler()
        console.setFormatter(color_formatter)

        self.addHandler(console)
        return
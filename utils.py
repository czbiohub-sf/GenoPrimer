import primer3
import requests
from Bio.Seq import Seq
import re
import csv
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

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

def get_primers(inputSeq, left10kb, right10kb, prod_size_lower, prod_size_upper, num_return, step_size):
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
            'SEQUENCE_INCLUDED_REGION': [1 + step_size*3, (len(inputSeq)-1) - step_size*3] # 90bp buffer zone on each side
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
    User_dict1["SEQUENCE_INCLUDED_REGION"] = [User_dict1["SEQUENCE_INCLUDED_REGION"][0] - step_size,
                                              User_dict1["SEQUENCE_INCLUDED_REGION"][1] + step_size]
    User_dict2["PRIMER_PRODUCT_SIZE_RANGE"] = [[User_dict2["PRIMER_PRODUCT_SIZE_RANGE"][0][0] + step_size*2,
                                                 User_dict2["PRIMER_PRODUCT_SIZE_RANGE"][0][1] + step_size*2]]

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
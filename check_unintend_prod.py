from BLAST_utils import check_blastDB_human
from subprocess import Popen
import os
import pandas as pd

def check_unintended_products(dict_primers, len_input, cut_chr , cut_coord, nonspecific_primers,fhlog):
    """
    checks unintended_products and remove primers having unintended products
    :param dict_primers:

    :param len_input: the length of primer
    :return: dict_primers(updated) #add a key PRIMER_PAIR_NUM_RETURNED_SPECIFIC
    """

    dict_primers["UNSPECIFIC_PRIMER_PAIR_idx"] = set()

    numThread2use = max([1, os.cpu_count()-2]) # for BLAST, use all CPUs except 2

    if dict_primers["PRIMER_PAIR_NUM_RETURNED"] == 0:
        return dict_primers

    max_pcr_prod_size = 6000
    print(f"Listing all possible PCR products < {max_pcr_prod_size} bp", flush=True)
    fhlog.write(f"Listing all possible PCR products < {max_pcr_prod_size} bp\n")

    dict_primer_len = {} #used to track the length of the primers, for 3' tracking purposes
    #create primer.fa
    empty_file_flag = 1 # this value =1 if primer.fa is empty
    tmp_fa = f"chr{cut_chr}_{cut_coord}_primer.fa"
    with open(tmp_fa, "w") as wfh:
        for i in range(0, dict_primers["PRIMER_PAIR_NUM_RETURNED"]):
            Lseq = dict_primers["PRIMER_LEFT_{}_SEQUENCE".format(i)]
            Rseq = dict_primers["PRIMER_RIGHT_{}_SEQUENCE".format(i)]
            #check if primers are in nonspecific list
            if not f"{Lseq}{Rseq}" in nonspecific_primers.keys():
                wfh.write(f">{i}_F\n{Lseq}\n>{i}_R\n{Rseq}\n")
                dict_primer_len[f"{i}_F"] = len(Lseq)
                dict_primer_len[f"{i}_R"] = len(Rseq)
                empty_file_flag = 0
            else:
                print(f"skipping {i}_F {i}_R due to nonspecific product detected in previous iteration(s)")
                fhlog.write(f"skipping {i}_F {i}_R due to nonspecific product detected in previous iteration(s)\n")
                dict_primers["UNSPECIFIC_PRIMER_PAIR_idx"].add(i)

    if empty_file_flag==0:
        #blast
        # check blastDB (human) and also return BLAST bin directory
        BLAST_bin, exe_suffix, BLAST_db_path = check_blastDB_human()
        # print(BLAST_bin)
        # print(exe_suffix)
        # print(BLAST_db_path)
        # specify BLAST db and query
        query = tmp_fa
        # start BLAST
        cmd = [f"{BLAST_bin}blastn{exe_suffix}", "-task", "blastn-short", "-query", f"{query}", "-db", f"{BLAST_db_path}", "-num_threads", f"{numThread2use}", "-perc_identity", "75", "-outfmt", "6 qseqid sseqid qstart qend sstart send pident mismatch", "-out", f"{query}.out"]
        p = Popen(cmd, universal_newlines=True)
        p.communicate()  # now wait for the process to finish
        os.remove(tmp_fa)

        #parse blast out

        df = pd.read_csv(f"{query}.out", sep = "\t", names=["qseqid", "sseqid", "qstart", "qend", "sstart", "send","pident", "mismatch"])

        # go through each primer pair
        primer_names = df["qseqid"].unique()
        primer_idx = set([i.split("_")[0] for i in df["qseqid"].unique()])
        primer_idx = sorted(primer_idx, key=float) #sort
        ## go through all primer pairs
        for idx in primer_idx:
            df_1pair = df[(df["qseqid"] == idx + "_F") | (df["qseqid"] == idx + "_R")]
            dict_primers = find_primer_parings(df_1pair=df_1pair,dict_primers=dict_primers, max_pcr_prod_size = max_pcr_prod_size, dict_primer_len=dict_primer_len, idx= idx, cut_chr =cut_chr, cut_coord = cut_coord, fhlog = fhlog)

        os.remove(f"{query}.out")
    return dict_primers


def find_primer_parings(df_1pair, dict_primers, max_pcr_prod_size, dict_primer_len, idx,cut_chr, cut_coord, fhlog):
    ## go through by each chr
    for current_chr in df_1pair["sseqid"].unique():
        df_1pair_1chr = df_1pair[df_1pair["sseqid"] == current_chr]
        #subset the dataframe (require primer 3' end to match), othewise the df may be very big
        F_len = dict_primer_len[f"{idx}_F"]
        R_len = dict_primer_len[f"{idx}_R"]
        df_1pair_1chr = df_1pair_1chr[((df_1pair_1chr["qseqid"] == f"{idx}_F") & (df_1pair_1chr["qend"] == F_len)) |
                                      ((df_1pair_1chr["qseqid"] == f"{idx}_R") & (df_1pair_1chr["qend"] == R_len))]
        current_pair_intended_prod_count = 0 #keep track of the number of intented targets for each primer pair

        ##go through all possible combinations that may produce PCR products
        combinations = []
        for idx_i, row_i in df_1pair_1chr.iterrows():
            for idx_j, row_j in df_1pair_1chr.iterrows():
                if not set([idx_i, idx_j]) in combinations:  # avoid i-j j-i repeats
                    combinations.append(set([idx_i, idx_j]))
                    strand_i = row_i["sstart"] - row_i["send"]
                    strand_j = row_j["sstart"] - row_j["send"]
                    coord_5p_i = row_i["sstart"]
                    coord_3p_i = row_i["send"]
                    coord_5p_j = row_j["sstart"]
                    coord_3p_j = row_j["send"]

                    # require matches be on different strands, and facing each other
                    if strand_i * strand_j < 0 and (
                            min([coord_5p_i, coord_5p_j]) < coord_3p_i < max([coord_5p_i, coord_5p_j])) and (
                            min([coord_5p_i, coord_5p_j]) < coord_3p_j < max([coord_5p_i, coord_5p_j])):
                        max_coord = max([row_i["sstart"], row_i["send"], row_j["sstart"], row_j["send"]])
                        min_coord = min([row_i["sstart"], row_i["send"], row_j["sstart"], row_j["send"]])
                        dist = max_coord - min_coord
                        # require the PCR product to be between 20bp and  max_pcr_prod_size
                        if 20 <= dist <= max_pcr_prod_size:
                            if float(row_i["pident"]) >= 80.0 and float(row_j["pident"]) >= 80.0:  # row_i["mismatch"]<=3 and row_j["mismatch"]<=3 and
                                # require less than 3 mismatches
                                if row_i["qend"] == dict_primer_len[row_i["qseqid"]] and row_j["qend"] == dict_primer_len[row_j["qseqid"]]:  # require the 3' end to match
                                    # check if this product is intended
                                    if str(current_chr) == str(cut_chr) and min_coord < cut_coord < max_coord:  # intended product
                                        # pass
                                        print(f"{row_i['qseqid']} chr={row_i['sseqid']} {row_i['sstart']}-{row_i['send']} {row_j['qseqid']}  chr={row_j['sseqid']} {row_j['sstart']}-{row_j['send']} product_size = {dist} (intended PCR product)",flush=True)
                                        fhlog.write(f"{row_i['qseqid']} chr={row_i['sseqid']} {row_i['sstart']}-{row_i['send']} {row_j['qseqid']}  chr={row_j['sseqid']} {row_j['sstart']}-{row_j['send']} product_size = {dist} (intended PCR product)\n")
                                        current_pair_intended_prod_count += 1
                                        if current_pair_intended_prod_count > 1: #flag primer pair unspecific if more than one intended PCR product is found
                                            dict_primers["UNSPECIFIC_PRIMER_PAIR_idx"].add(idx)
                                    else:  # flag primer for having unintened product
                                        # mark primers with unintended products
                                        dict_primers["UNSPECIFIC_PRIMER_PAIR_idx"].add(idx)
                                        print(f"{row_i['qseqid']} chr={row_i['sseqid']} {row_i['sstart']}-{row_i['send']} {row_j['qseqid']}  chr={row_j['sseqid']} {row_j['sstart']}-{row_j['send']} product_size = {dist} (*unintended* PCR product) Skipping listing other PCR products for this primer",flush=True)
                                        fhlog.write(f"{row_i['qseqid']} chr={row_i['sseqid']} {row_i['sstart']}-{row_i['send']} {row_j['qseqid']}  chr={row_j['sseqid']} {row_j['sstart']}-{row_j['send']} product_size = {dist} (*unintended* PCR product) Skipping listing other PCR products for this primer\n")
                                        return dict_primers # this ends the function, thus stops checking this primer any further
    return dict_primers



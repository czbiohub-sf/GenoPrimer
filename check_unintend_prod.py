from BLAST_utils import check_blastDB
from Bowtie_utils import check_BowtieIdx
from subprocess import Popen
from subprocess import PIPE
from subprocess import DEVNULL
import os
import pandas as pd

def check_unintended_products(dict_primers, len_input, ref, cut_chr , cut_coord, nonspecific_primers,thread, fhlog,outdir, aligner="Bowtie"):
    """
    checks unintended_products and remove primers having unintended products
    :param dict_primers:
    :param Aligner: BLAST or Bowtie

    :param len_input: the length of primer
    :return: dict_primers(updated) #add a key PRIMER_PAIR_NUM_RETURNED_SPECIFIC
    """

    dict_primers["UNSPECIFIC_PRIMER_PAIR_idx"] = set()

    if thread == "auto":
        numThread2use = max([1, os.cpu_count()-2]) # for BLAST, use all CPUs except 2
    else:
        numThread2use = int(thread)

    if dict_primers["PRIMER_PAIR_NUM_RETURNED"] == 0:
        return dict_primers

    max_pcr_prod_size = 6000
    max_genomic_hits = 5000 #max number of hsps across all hits allowed for each primer, offending primers will not be considered

    #print(f"Listing all possible PCR products < {max_pcr_prod_size} bp", flush=True)
    fhlog.write(f"Listing all possible PCR products < {max_pcr_prod_size} bp\n")

    dict_primer_len = {} #used to track the length of the primers, for 3' tracking purposes
    #create primer.fa
    empty_file_flag = 1 # this value =1 if primer.fa is empty
    tmp_fa = f"chr{cut_chr}_{cut_coord}_primer.fa"
    tmp_fa = os.path.join(outdir,tmp_fa)
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
                #print(f"skipping {i}_F {i}_R due to nonspecific product detected in previous iteration(s)")
                fhlog.write(f"skipping {i}_F {i}_R due to nonspecific product detected in previous iteration(s)\n")
                dict_primers["UNSPECIFIC_PRIMER_PAIR_idx"].add(i)

    if empty_file_flag==0:

        if aligner == "BLAST" :
            #blast
            # check blastDB (human) and also return BLAST bin directory
            BLAST_bin, exe_suffix, BLAST_db_path = check_blastDB(ref)
            # print(BLAST_bin)
            # print(exe_suffix)
            # print(BLAST_db_path)
            # specify BLAST db and query
            query = tmp_fa
            # start BLAST
            cmd = [f"{BLAST_bin}blastn{exe_suffix}", "-task", "blastn-short", "-query", f"{query}", "-db", f"{BLAST_db_path}", "-max_hsps", "2000", "-num_threads", f"{numThread2use}", "-perc_identity", "75", "-outfmt", "6 qseqid sseqid qstart qend sstart send pident mismatch", "-out", f"{query}.out"]
            p = Popen(cmd, universal_newlines=True)
            p.communicate()  # now wait for the process to finish
            os.remove(tmp_fa)

            #pre-parse alignment out, flag primers with too many hits in the genome
            uniq_cout = f"{query}.out.uniq_count"
            with open(uniq_cout, "w") as f, open(f"{query}.out", "r") as i:
                uniq_count = dict()
                for line in i:
                    pname = line.split("\t")[0]
                    if pname in uniq_count.keys():
                        uniq_count[pname] += 1
                    else:
                        uniq_count[pname] = 1
                for k in uniq_count.keys():
                    f.write(f"{uniq_count[k]}\t{k}\n")

            to_remove = []
            with open(uniq_cout, "r") as f:
                for line in f:
                    line = line.lstrip().rstrip()
                    num,primer = line.split("\t")
                    #print(f"{num},{primer}")
                    if int(num) >= max_genomic_hits:
                        to_remove.append(primer) #remove primers with too many hits in the genome
            #print(f"to_remove:{to_remove}")

            os.remove(uniq_cout)

            with open(f"{query}.out", "r") as f, open(f"{query}.out.preparsed", "w") as w:
                for line in f:
                    primer = line.split("\t")[0]
                    if not primer in to_remove:
                        w.write(line)

            os.remove(f"{query}.out")

        elif aligner == "Bowtie":
            Bowtie_bin, exe_suffix, Bowtie_idx_path = check_BowtieIdx(ref)
            query = tmp_fa
            # run Bowtie
            with open(f"{query}.out", "w") as outfile:
                cmd = [f"{Bowtie_bin}bowtie{exe_suffix}", "-f", "-k", "1000", "-v", "3", "--threads", "4", "--suppress", "5,6", "--quiet", "-x", f"{Bowtie_idx_path}", f"{query}"]
                p = Popen(cmd, universal_newlines=True, stdout=outfile, stderr=DEVNULL)
                p.communicate()  # now wait for the process to finish
                os.remove(tmp_fa)

            # parse Bowtie out, flag primers with too many hits in the genome, remove alignments that have mismatches after the 15th base
            to_remove = []
            with open(f"{query}.out.preparsed", "w") as f, open(f"{query}.out", "r") as r:
                uniq_count = dict()
                for line in r:
                    fields = line.strip().split("\t")
                    pname = fields[0]
                    largest_mismatch_pos = 0
                    num_missmatch = 0
                    #check mismatch_desc fields to update largest_mismatch_pos
                    if len(fields) > 5:
                        mismatch_desc = fields[5]
                        mismatch_desc_fields = mismatch_desc.split(",")           
                        num_missmatch = len(mismatch_desc_fields)
                        if len(mismatch_desc_fields) > 0: #mismatch_desc string is not empty
                            for c in mismatch_desc_fields:
                                mismatch_pos = int(c.split(":")[0])
                                if mismatch_pos > largest_mismatch_pos:
                                    largest_mismatch_pos = mismatch_pos
                    if largest_mismatch_pos <= 15: # require no mismatch after the 15th base
                        if int(fields[4]) <= max_genomic_hits: # require no more than max_genomic_hits hits
                            strand = fields[1]
                            if strand == "+":
                                sstart = int(fields[3])
                                send = int(fields[3]) + 20 # use avg primer length of 20
                            else:
                                sstart = int(fields[3]) + 20 # use avg primer length of 20
                                send = int(fields[3])
                            pident = 100 - num_missmatch/20*100
                            f.write(f"{pname}\t{fields[2]}\t1\t20\t{sstart}\t{send}\t{pident}\t{num_missmatch}\n")
                        else:
                            to_remove.append(pname) #remove primers with too many hits in the genome
                os.remove(f"{query}.out")

        else:
            sys.exit(f"[ERROR] Invalid choice of aligner: {aligner}")

        #flag all the removed-primers in dict_primers as UNSPECIFIC_PRIMER_PAIR
        for primer in to_remove:
            idx = primer.split("_")[0]
            dict_primers["UNSPECIFIC_PRIMER_PAIR_idx"].add(idx)
            #print(f"idx:{idx} removed")
        
        #parse blast out
        df = pd.read_csv(f"{query}.out.preparsed", sep = "\t", names=["qseqid", "sseqid", "qstart", "qend", "sstart", "send","pident", "mismatch"])
        os.remove(f"{query}.out.preparsed")

        # go through each primer pair
        primer_names = df["qseqid"].unique()
        primer_idx = set([i.split("_")[0] for i in df["qseqid"].unique()])
        primer_idx = sorted(primer_idx, key=float) #sort
        ## go through all primer pairs
        for idx in primer_idx:
            df_1pair = df[(df["qseqid"] == idx + "_F") | (df["qseqid"] == idx + "_R")]
            dict_primers = find_primer_parings(df_1pair=df_1pair,dict_primers=dict_primers, max_pcr_prod_size = max_pcr_prod_size, dict_primer_len=dict_primer_len, idx= idx, cut_chr =cut_chr, cut_coord = cut_coord, fhlog = fhlog)
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
            #get subset of primers that are within [max_pcr_prod_size] distance away from the current primer
            anchor = row_i["sstart"]
            df_1pair_1chr_subset = df_1pair_1chr[(df_1pair_1chr["sstart"] >= (anchor - max_pcr_prod_size - 100)) & (df_1pair_1chr["sstart"] <= (anchor + max_pcr_prod_size + 100))]
            for idx_j, row_j in df_1pair_1chr_subset.iterrows():
                if not set([idx_i, idx_j]) in combinations:  # avoid i-j j-i repeats
                    combinations.append(set([idx_i, idx_j]))
                    strand_i = row_i["sstart"] - row_i["send"]
                    strand_j = row_j["sstart"] - row_j["send"]
                    coord_5p_i = row_i["sstart"]
                    coord_3p_i = row_i["send"]
                    coord_5p_j = row_j["sstart"]
                    coord_3p_j = row_j["send"]

        # for idx_i, row_i in df_1pair_1chr.iterrows():
        #     for idx_j, row_j in df_1pair_1chr.iterrows():
        #         if not set([idx_i, idx_j]) in combinations:  # avoid i-j j-i repeats
        #             combinations.append(set([idx_i, idx_j]))
        #             strand_i = row_i["sstart"] - row_i["send"]
        #             strand_j = row_j["sstart"] - row_j["send"]
        #             coord_5p_i = row_i["sstart"]
        #             coord_3p_i = row_i["send"]
        #             coord_5p_j = row_j["sstart"]
        #             coord_3p_j = row_j["send"]

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
                                        #print(f"{row_i['qseqid']} chr={row_i['sseqid']} {row_i['sstart']}-{row_i['send']} {row_j['qseqid']}  chr={row_j['sseqid']} {row_j['sstart']}-{row_j['send']} product_size = {dist} (intended PCR product)",flush=True)
                                        fhlog.write(f"{row_i['qseqid']} chr={row_i['sseqid']} {row_i['sstart']}-{row_i['send']} {row_j['qseqid']}  chr={row_j['sseqid']} {row_j['sstart']}-{row_j['send']} product_size = {dist} (intended PCR product)\n")
                                        current_pair_intended_prod_count += 1
                                        if current_pair_intended_prod_count > 1: #flag primer pair unspecific if more than one intended PCR product is found
                                            dict_primers["UNSPECIFIC_PRIMER_PAIR_idx"].add(idx)
                                    else:  # flag primer for having unintened product
                                        # mark primers with unintended products
                                        dict_primers["UNSPECIFIC_PRIMER_PAIR_idx"].add(idx)
                                        #print(f"{row_i['qseqid']} chr={row_i['sseqid']} {row_i['sstart']}-{row_i['send']} {row_j['qseqid']}  chr={row_j['sseqid']} {row_j['sstart']}-{row_j['send']} product_size = {dist} (*unintended* PCR product) Skipping listing other PCR products for this primer",flush=True)
                                        fhlog.write(f"{row_i['qseqid']} chr={row_i['sseqid']} {row_i['sstart']}-{row_i['send']} {row_j['qseqid']}  chr={row_j['sseqid']} {row_j['sstart']}-{row_j['send']} product_size = {dist} (*unintended* PCR product) Skipping listing other PCR products for this primer\n")
                                        return dict_primers # this ends the function, thus stops checking this primer any further
    return dict_primers



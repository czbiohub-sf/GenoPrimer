import pandas as pd
import os, sys, pickle

input_file = "input/OC1_protospacer.coord.multi-resolved.csv"
ENST_info = "BLAST_databases/ENST_info.pickle"

def read_pickle_files(file):
    if os.path.isfile(file):
        with open(file, "rb") as handle:
            mydict = pickle.load(handle)
        return mydict
    else:
        sys.exit(f"Cannot open file: {file}")

def load_ENST_info():   
    #load ENST_info
    print("loading ENST gene model info")
    ENST_info = read_pickle_files(os.path.join("BLAST_databases",'ENST_info.pickle'))
    return ENST_info

def check_ATG_at_exonEnd(my_transcript):
    """
    input: transcript object
    output: Bool
    """
    transcript_type = my_transcript.description.split("|")[1]
    if transcript_type == "protein_coding":  # only look at protein-coding transcripts
        # constructing the list of cds
        cdsList = [feat for feat in my_transcript.features if feat.type == "CDS"]
        if len(cdsList) >= 1:  # has more than 1 cds
            CDS_first = cdsList[0]
            cds_len = (
                abs(CDS_first.location.start - CDS_first.location.end) + 1
            )  # for ATG to be the end of the exon, the first exon length is 3bp
            if cds_len == 3:
                return True
            else:
                return False
        else:
            return False
    else:
        return False
    
def get_start_stop_loc(ENST_ID, ENST_info):
    """
    Get the chromosomal location of start and stop codons
    input: ENST_ID, ENST_info
    output: a list of three items
            ATG_loc: [chr,start,end,strand]  #start < end
            stop_loc: [chr,start,end,strand] #start < end
            Exon_end_ATG: Bool
    """
    my_transcript = ENST_info[ENST_ID]  # get the seq record
    # constructing the list of cds
    cdsList = [feat for feat in my_transcript.features if feat.type == "CDS"]
    CDS_first = cdsList[0]
    CDS_last = cdsList[len(cdsList) - 1]
    # check if ATG is at the end of the first exon
    if check_ATG_at_exonEnd(my_transcript):
        CDS_first = cdsList[
            1
        ]  # use the second cds if ATG is at the end of the first exon
    # get start codon location
    if CDS_first.strand == 1:
        ATG_loc = [
            CDS_first.location.ref,
            CDS_first.location.start + 0,
            CDS_first.location.start + 2,
            1,
        ]  # format [start, end, strand]
    else:
        stop_loc = [
            CDS_first.location.ref,
            CDS_first.location.start + 0,
            CDS_first.location.start + 2,
            -1,
        ]
    # get stop codon location
    if CDS_last.strand == 1:
        stop_loc = [
            CDS_last.location.ref,
            CDS_last.location.end - 2,
            CDS_last.location.end + 0,
            1,
        ]
    else:
        ATG_loc = [
            CDS_last.location.ref,
            CDS_last.location.end - 2,
            CDS_last.location.end + 0,
            -1,
        ]

    return [ATG_loc, stop_loc]


def get_end_pos_of_ATG(ATG_loc):
    """
    input: ATG_loc              [chr,start,end,strand] #start < end
    output: the pos of G in ATG [chr,pos,strand]       #start < end
    """
    strand = ATG_loc[3]
    if str(strand) == "+" or str(strand) == "1":
        return [ATG_loc[0], ATG_loc[2], ATG_loc[3]]
    else:
        return [ATG_loc[0], ATG_loc[1]-1, ATG_loc[3]]  # -1 because we want the insertion site to be downstream of the pos
    
def get_start_pos_of_stop(stop_loc):
    """
    input: stop_loc                                 [chr,start,end,strand]  #start < end
    output: the pos of first base in the stop codon [chr,pos,strand]        #start < end
    """
    strand = stop_loc[3]
    if str(strand) == "+" or str(strand) == "1":
        return [stop_loc[0], stop_loc[1]-1, stop_loc[3]] # -1 because we want the insertion site to be downstream of the pos
    else:
        return [stop_loc[0], stop_loc[2], stop_loc[3]]
    

if __name__ == "__main__":
    ENST_info = load_ENST_info()
    df = pd.read_csv(input_file)
    with open(os.path.join(f"{input_file}.cut2insert.csv"), 'w') as outcsv:
        outcsv.write(",".join(list(df.columns) + ["insert_position","cut-to-insert distance"])) #header
        outcsv.write("\n")

        for index, row in df.iterrows():
            ENST_ID = row["ENST"]
            print(ENST_ID)
            csvrow = [str(item) for item in row]
            if ENST_ID in ENST_info.keys():
                ATG_loc, stop_loc = get_start_stop_loc(ENST_ID, ENST_info)
                terminus = row["terminus_to_tag"].upper()
                if terminus == "N":
                    end_of_ATG_loc = get_end_pos_of_ATG(ATG_loc)
                    insert_pos = end_of_ATG_loc[1]
                    if end_of_ATG_loc[0] != row["Ensemble_chr"]:
                        print(f"chr doesn't match for insert site and cutsite")
                        outcsv.write(",".join(csvrow) + ",chr doesn't match for insert site and cutsite\n") 
                        continue
                elif terminus == "C":
                    start_of_stop_loc = get_start_pos_of_stop(stop_loc)
                    insert_pos = start_of_stop_loc[1]
                    if start_of_stop_loc[0] != row["Ensemble_chr"]:
                        print(f"chr doesn't match for insert site and cutsite")
                        outcsv.write(",".join(csvrow) + ",chr doesn't match for insert site and cutsite\n") 
                        continue
                elif terminus == "INT":
                    outcsv.write(",".join(csvrow) + ",terminues_to_tag is INT\n") 
                    continue
                else:
                    sys.exit(f"unknown terminus {terminus} terminus_to_tag must be either N or C")
                
                cut_to_insert_distance = row["gRNACut_in_chr"] - insert_pos
                outcsv.write(",".join(csvrow) + f",{insert_pos},{cut_to_insert_distance}\n")
            else:
                #ENST not found
                outcsv.write(",".join(csvrow) + ",Cannot find ENST ID in the database\n") 


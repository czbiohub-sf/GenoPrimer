import argparse
import sys
import linecache
import primer3
import csv
import pandas as pd

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

#####################
##      main       ##
#####################    
def main():
    try: 

            df = pd.read_excel(os.path.join(config['csv']), sheet_name="Oocyst_count")

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

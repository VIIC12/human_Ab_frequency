#!/usr/bin/env python
import argparse
import logging

def error_checker(args):
    """
    Checks user input for typical incorrect search queries (like forbidden characters) and prints out error message if necessary.
    """
    import os
    import re
    string = re.compile('[@_!#$%^&()<>-?/\}{~:a-zA-Z]')
    seq = re.compile('[@_!#$%^&()<>-?/\|}{~:0-9BJOUXZ]')
    if args.overwrite == 0:
        try:
            os.makedirs(args.outputdir)
        except FileExistsError:
            logging.error('Output directory already exists & overwrite is disabled.')
            return False

    if string.search(args.IGHV or args.IGHD or args.IGHJ) != None:
        logging.error("User input error: Do only enter numbers, - and * for IGHV/D/J!")
        return False
    
    if seq.search(args.h3_motif) != None:
        logging.error("User input error: Do only enter numbers, - and * for h3_motif!")
        return False
    
    if args.CDRH3_length != None:
        if len(args.h3_motif) > int(args.CDRH3_length):
            logging.error("User input error: The CDRH3 sequence length is longer than the CDRH3 length!")
            return False
    
    return True

def adjust_gene(gene):
    adjusted_parts = []
    parts = gene.split("|")
    for part in parts:
        if "-" in part:
            adjusted_parts.append(part)
        elif "-" not in part:
            adjusted_parts.append(part + "-")
    return "|".join(adjusted_parts)

def adjust_input(args):
    """
    Adjust user input
    """
    if args.IGHV:
        args.IGHV = adjust_gene(args.IGHV)
    if args.IGHD:
        args.IGHD = adjust_gene(args.IGHD)
    if args.IGHJ:
        adjusted_parts = []
        parts = args.IGHJ.split("|")
        for part in parts:
            adjusted_parts.append("J" + part)
        args.IGHJ = "|".join(adjusted_parts)
    if args.CDRH3_length is not None:
        args.CDRH3_length = str(args.CDRH3_length)
    else:
        args.CDRH3_length = ""

    if args.database[-1] != '/':
        args.database = args.database + '/'
    if args.outputdir[-1] != '/':
        args.outputdir = args.outputdir + '/'
    return args

def main(args):
    """
    Human Ab freuqency app
    """
    import pandas as pd
    import glob
    import json
    import gzip
    import numpy as np
    import time
    from progress.bar import Bar

    logging.basicConfig(level=logging.INFO, format='%(message)s')
    logging.getLogger('numexpr').setLevel(logging.WARNING)
    start_time = time.time()

    logging.info(
        """
        ╔─────────────────────────────────────────────────────────────────────────────╗
        │██╗  ██╗██╗   ██╗███╗   ███╗ █████╗ ███╗   ██╗     █████╗ ██████╗            │
        │██║  ██║██║   ██║████╗ ████║██╔══██╗████╗  ██║    ██╔══██╗██╔══██╗           │
        │███████║██║   ██║██╔████╔██║███████║██╔██╗ ██║    ███████║██████╔╝           │
        │██╔══██║██║   ██║██║╚██╔╝██║██╔══██║██║╚██╗██║    ██╔══██║██╔══██╗           │
        │██║  ██║╚██████╔╝██║ ╚═╝ ██║██║  ██║██║ ╚████║    ██║  ██║██████╔╝           │
        │╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚═╝╚═╝  ╚═╝╚═╝  ╚═══╝    ╚═╝  ╚═╝╚═════╝            │
        │                                                                             │
        │███████╗██████╗ ███████╗ ██████╗ ██╗   ██╗███████╗███╗   ██╗ ██████╗██╗   ██╗│
        │██╔════╝██╔══██╗██╔════╝██╔═══██╗██║   ██║██╔════╝████╗  ██║██╔════╝╚██╗ ██╔╝│
        │█████╗  ██████╔╝█████╗  ██║   ██║██║   ██║█████╗  ██╔██╗ ██║██║      ╚████╔╝ │
        │██╔══╝  ██╔══██╗██╔══╝  ██║▄▄ ██║██║   ██║██╔══╝  ██║╚██╗██║██║       ╚██╔╝  │
        │██║     ██║  ██║███████╗╚██████╔╝╚██████╔╝███████╗██║ ╚████║╚██████╗   ██║   │
        │╚═╝     ╚═╝  ╚═╝╚══════╝ ╚══▀▀═╝  ╚═════╝ ╚══════╝╚═╝  ╚═══╝ ╚═════╝   ╚═╝   │
        │                                                                             │
        │ █████╗ ██████╗ ██████╗                                                      │
        │██╔══██╗██╔══██╗██╔══██╗                                                     │
        │███████║██████╔╝██████╔╝                                                     │
        │██╔══██║██╔═══╝ ██╔═══╝                                                      │
        │██║  ██║██║     ██║      Version 1.0                                         │
        │╚═╝  ╚═╝╚═╝     ╚═╝      ©2023 T.U. Schlegel (tom.schlegel@uni-leipzig.de)   │
        ╚─────────────────────────────────────────────────────────────────────────────╝
        """
    )

    #? Abort if input errors
    if not error_checker(args):
        return
    logging.debug("User input seems to be correct.")    

    #? Adjust User input
    args = adjust_input(args)
    logging.info(f"Your search request:\n \
        IGHV: {args.IGHV}\n \
        IGHD: {args.IGHD}\n \
        IGHJ: {args.IGHJ}\n \
        CDRH3_length: {args.CDRH3_length}\n \
        h3_motif: {args.h3_motif}\n \
        database: {args.database}\n \
        outputdir: {args.outputdir}\n \
        full_results: {args.full_results}\n \
        overwrite: {args.overwrite}\n")
    all_trial_info = dict()
    dffull_sort = pd.DataFrame()

    db = args.database
    with Bar('Processing...', max = len([file for file in glob.glob(db+'*.csv.gz')])) as bar:
        for file in glob.glob(db+'*.csv.gz'):
            logging.debug(f"Reading file: {file}")  
            #df_single = pd.read_csv(file, skiprows=1)
            df_single = pd.read_csv(file, skiprows=1, usecols = ['v_call','d_call','j_call','junction_aa_length','junction_aa'])
            dfsingle_sort = df_single.loc[
            df_single['v_call'].astype(str).str.contains(args.IGHV) & \
            df_single['d_call'].astype(str).str.contains(args.IGHD) & \
            df_single['j_call'].astype(str).str.contains(args.IGHJ) & \
            df_single['junction_aa_length'].astype(str).str.startswith(args.CDRH3_length) & \
            df_single['junction_aa'].astype(str).str.contains(args.h3_motif)]

            rowstotal = len(dfsingle_sort.index)
            value = rowstotal

            #? Add Meta-info
            with gzip.open(file, 'rt', encoding='utf-8') as file:
                trial_info = json.loads(file.readline().replace('"{', '{').replace('}"', '}').replace('""', '"'))

            sumtotal = int(trial_info['Unique sequences'])
                    
            #? Update dic
            # if Subject already as key in dict, than add new value and total to exisiting value        
            if trial_info['Subject'] in all_trial_info:
                newvalue = all_trial_info.get(trial_info['Subject'])[0] + value
                newtotal = all_trial_info.get(trial_info['Subject'])[1] + sumtotal
                all_trial_info.update({trial_info['Subject']: [newvalue, newtotal]})
                logging.debug("Subject already in dict - add hits and total")  
            else:
                all_trial_info[trial_info['Subject']] = [value, sumtotal]
                logging.debug("Subject not in dict - create new Subject with hits and total")  

            dfsingle_sort = dfsingle_sort.reset_index(drop=True)
            dfsingle_sort['Subject source'] = trial_info['Subject']
            dffull_sort = pd.concat([dffull_sort, dfsingle_sort])
            bar.next()

    resdf = pd.DataFrame.from_dict(all_trial_info).rename(index={0 :'Hits', 1:'Total sequences'})
    resdf.loc['Percentage of hits'] = (resdf.iloc[0] / resdf.iloc[1] *100).round(2)
    resdf.loc['Hits per 1 million'] = (1000000*resdf.iloc[0] / resdf.iloc[1]).round(0)
    resdf['Meta'] = np.nan
    resdf.iloc[0, resdf.columns.get_loc('Meta')] = f"IGHV: {args.IGHV}, IGHD: {args.IGHD}, IGHJ: {args.IGHJ}, CDRH3-length: {args.CDRH3_length}, H3-motif: {args.h3_motif}"
    resdf.to_csv(args.outputdir+"results.csv")

    if args.full_results == 1:
        logging.info("Merging fullresults file... (This may take a while.)")
        dffull_sort.to_csv(args.outputdir+"fullresults.csv")
    logging.info(f"Done! Search duration: {time.time() - start_time:.2f} sec. Thank you for using the Human Ab Frequency App!")

if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--IGHV", type=str, default="", help='"3-" or "3-22" or multiple by using | like: "3-20|3-22"')
    argparser.add_argument("--IGHD", type=str, default="", help='"3-" or "3-22" or multiple by using | like: "3-20|3-22"')
    argparser.add_argument("--IGHJ", type=str, default="", help='"5" or multiple by using | like: "4|5"')
    argparser.add_argument("--CDRH3_length", type=int, default=None, help='Length of CDHR3 region as int (WITH C-X-W, so if necessary add 2)')
    argparser.add_argument("--h3_motif", type=str, default="", help='"." for one, ".*" for 0-many, like "YY.D.*G"')
    argparser.add_argument("--database", type=str, default="", help='Set database directory')   
    argparser.add_argument("--full_results", type=int, default="1", help='Safe full results table. 1 for True, 0 for False')
    argparser.add_argument("--outputdir", type=str, default="output/", help='Output directory')   
    argparser.add_argument("--overwrite", type=int, default="0", help='1 for True, 0 for False')   
    args = argparser.parse_args() 
    main(args)

__author__ = "Tom U. Schlegel"
__contact__ = "tom.schlegel@uni-leipzig.de"
__license__ = "GNU GPLv3"
__version__ = "1.0"
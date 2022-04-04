#!/usr/bin/env python

import pandas as pd
import glob
import argparse
import os

def arg_parser():
    """Handles the argument in and output on the command line, returns the
    arguments given by the user"""
    argp = argparse.ArgumentParser(description="Aggregates all UMI_counttables to a single overview")
    argp.add_argument("-i", "--inputDir", default=str,
                      help="""location containing all UMI_counttable.csv files to parse""")
    argp.add_argument("-o", "--outputDir", default=str, nargs="?",
                      help="""output file of the program, default = stdout""")
    return argp.parse_args()

def summary(ls, column, outputDir):

    frames = []

    # loop and open in panda
    for f in ls:
        print(f)
        name = os.path.basename(f)
        try:
            df = pd.read_csv(f, sep=',', engine='python', comment='##',
                             header=0)
            print(df)
            df["amplicon"] = df["fw_name"] + "|" + df["rv_name"]
            df = df[["amplicon",column]]
            df = df.groupby(['amplicon']).sum()
            print(df)
            df = df.rename(columns={column: name})
            print(df)
        except pd.errors.EmptyDataError:
            continue
        frames.append(df)

    df_final = pd.concat(frames, axis=1)
    df_final = df_final.fillna(0)

    # store result
    df_final.to_csv(f"{outputDir}/UMI_counttable_total_{column}.txt", sep="\t")

def main():

    args = arg_parser()

    ls = glob.glob(f"{args.inputDir}/*_UMI_counttable.csv")
    print(ls)

    summary(ls, "percentage", args.outputDir)
    summary(ls, "Count", args.outputDir)

if __name__ == "__main__":
    # load arguments set global workdir
    # args = parse_args()
    # fill_html(args)
    print("Start")
    main()
    print("Finished")
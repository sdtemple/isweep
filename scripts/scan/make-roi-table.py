# Append info to region of interest.
# Seth Temple, GitHub: sdtemple
# Last modified: 2024-07-22

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import floor, ceil
import sys

def main():
    # Set up argparse
    parser = argparse.ArgumentParser(
        description="Filter, append, format info to region of interest."
    )

    parser.add_argument(
        '--input_excess_file',
        type=str,
        required=True,
        help="Input file containing excess IBD regions data."
    )

    parser.add_argument(
        '--input_scan_file',
        type=str,
        required=True,
        help="Prefix of chromosome files."
    )

    parser.add_argument(
        '--output_file',
        type=str,
        required=True,
        help="Output file for saving region of interest."
    )

    # parser.add_argument(
    #     '--input_prefix',
    #     type=str,
    #     required=True,
    #     help="Prefix of chromosome files."
    # )

    # parser.add_argument(
    #     '--input_suffix',
    #     type=str,
    #     required=True,
    #     help="Suffix of chromosome files."
    # )
    
    parser.add_argument(
        '--cM_cover', 
        type=float,
        default=1.0, 
        help="(default: 1.0) Threshold for covering centiMorgans (cM)."
    )
    
    parser.add_argument(
        '--cM_small', 
        type=float,
        default=2.0, 
        help="(default: 2.0) Threshold for small centiMorgans (cM)."
    )
    
    parser.add_argument(
        '--Mb_buffer', 
        type=float,
        default=2.0, 
        help="(default: 2) Buffer size in megabases (Mb)."
    )

    parser.add_argument(
        '--statistic', 
        type=str,
        default='COUNT', 
        help="(default: COUNT) Test statistic that exceeds threshold."
    )

    parser.add_argument(
        '--sweep', 
        type=int,
        default=1, 
        help="(default: 1) If > 0 (selection scan), format columns with model and alpha. If 0 (case scan), do not have these columns."
    )
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Extracting arguments
    cmcover = args.cM_cover
    cmsmall = args.cM_small
    mbbuffer = args.Mb_buffer * 1000000
    
    filein = f"{args.input_excess_file}"
    fileout = f"{args.output_file}"

    statistic = args.statistic
    minstat = 'MIN' + statistic
    maxstat = 'MAX' + statistic
    
    roitab = pd.read_csv(filein, sep='\t')
    bpcenter = []
    cmcenter = []

    # if there is no excess IBD make a blank file
    if roitab.shape[0] == 0:
        with open(fileout, 'w') as f:
            f.write('NAME\tCHROM\t'+ minstat +'\t' + maxstat + '\tBPCENTER\tCMCENTER\tBPLEFTCENTER\tBPRIGHTCENTER\n')
        sys.exit()

    # make magic happen

    filein1 = args.input_scan_file
    tab = pd.read_csv(filein1, sep='\t')

    # formatting
    for i in range(roitab.shape[0]):
        # reading
        rowing = roitab.iloc[i]
        left = floor(float(rowing.BPLEFT))
        right = ceil(float(rowing.BPRIGHT))
        chrom = floor(float(rowing.CHROM))
        subtab = tab[(tab['BPWINDOW'] >= left) & 
                     (tab['BPWINDOW'] <= right) & 
                     (tab['CHROM'] == chrom)]
        subtab['WEIGHT'] = subtab[statistic] / subtab[statistic].sum()
        subtab['CUMSUM'] = np.cumsum(subtab['WEIGHT'])

        # central tendency
        shape0 = subtab.shape[0]
        mx = 0
        moCM = 0
        moBP = 0
        for j in range(shape0):
            row = subtab.iloc[j]
            if row[statistic] > mx:
                moCM = row['CMWINDOW']
                moBP = row['BPWINDOW']
                mx = row[statistic]
        moBP = int(moBP)
        meCM = (subtab['WEIGHT'] * subtab['CMWINDOW']).sum()
        meBP = (subtab['WEIGHT'] * subtab['BPWINDOW']).sum()
        meBP = int(meBP)
        mdCM = subtab[subtab['CUMSUM'] >= 0.5]['CMWINDOW'].tolist()[0]
        mdBP = subtab[subtab['CUMSUM'] >= 0.5]['BPWINDOW'].tolist()[0]
        bpcenter.append(moBP)
        cmcenter.append(moCM)

    # adding bp, cm centrality to roi table
    roitab['BPCENTER'] = bpcenter
    roitab['CMCENTER'] = cmcenter

    roitab = roitab[(roitab['CMRIGHT'] - roitab['CMLEFT']) >= cmcover]
    roitab['BPLEFTCENTER'] = roitab['BPLEFT']
    roitab['BPRIGHTCENTER'] = roitab['BPRIGHT']
    roitab.loc[(roitab['CMRIGHT'] - roitab['CMLEFT']) <= cmsmall, 'BPLEFTCENTER'] = roitab.loc[(roitab['CMRIGHT'] - roitab['CMLEFT']) <= cmsmall, 'BPLEFTCENTER'] - mbbuffer
    roitab.loc[(roitab['CMRIGHT'] - roitab['CMLEFT']) <= cmsmall, 'BPRIGHTCENTER'] = roitab.loc[(roitab['CMRIGHT'] - roitab['CMLEFT']) <= cmsmall, 'BPRIGHTCENTER'] + mbbuffer
    roitab['BPLEFTCENTER'] = roitab['BPLEFTCENTER'].clip(lower=1)

    # sorting, giving generic names
    initcol = list(roitab.columns)
    roitab.sort_values(by=maxstat, ascending=False, inplace=True)
    nrow = roitab.shape[0]
    roitab['NAME'] = [f'hit{i}' for i in range(1, nrow+1)]
    if args.sweep > 0:
        finacol = ['NAME'] + initcol + ['MODEL', 'ALPHA']
        roitab['MODEL'] = ['a' for _ in range(1, nrow+1)]
        roitab['ALPHA'] = [0.95 for _ in range(1, nrow+1)]
    elif args.sweep <= 0:
        finacol = ['NAME'] + initcol
    roitab = roitab[finacol]
    if args.sweep > 0:
        roitab.to_csv(fileout, sep='\t', index=False, 
                      columns=['NAME', 'CHROM', minstat, maxstat,
                                'BPCENTER', 'CMCENTER', 'BPLEFTCENTER', 'BPRIGHTCENTER',
                                'MODEL', 'ALPHA'])
    elif args.sweep <= 0:
        roitab.to_csv(fileout, sep='\t',index=False,
                      columns=['NAME','CHROM', minstat, maxstat,
                               'BPCENTER','CMCENTER','BPLEFTCENTER','BPRIGHTCENTER'
                               ]
                      )

if __name__ == "__main__":
    main()
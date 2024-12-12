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
        '--input_file',
        type=str,
        required=True,
        help="Input file containing IBD regions data."
    )

    parser.add_argument(
        '--output_file',
        type=str,
        required=True,
        help="Output file for saving region of interest."
    )

    parser.add_argument(
        '--input_prefix',
        type=str,
        required=True,
        help="Prefix of chromosome files."
    )

    parser.add_argument(
        '--input_suffix',
        type=str,
        required=True,
        help="Suffix of chromosome files."
    )
    
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
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Extracting arguments
    cmcover = args.cM_cover
    cmsmall = args.cM_small
    mbbuffer = args.Mb_buffer * 1000000
    
    filein = f"{args.input_file}"
    filepre = f"{args.input_prefix}"
    filesuf = f"{args.input_suffix}"
    fileout = f"{args.output_file}"
    
    roitab = pd.read_csv(filein, sep='\t')
    bpcenter = []
    cmcenter = []

    # if there is no excess IBD make a blank file
    if roitab.shape[0] == 0:
        with open(fileout, 'w') as f:
            f.write('NAME\tCHROM\tMINIBD\tMAXIBD\tBPCENTER\tCMCENTER\tBPLEFTCENTER\tBPRIGHTCENTER\n')
        sys.exit()

    # make magic happen

    # formatting
    for i in range(roitab.shape[0]):
        # reading
        rowing = roitab.iloc[i]
        filein1 = f"{filepre}{int(float(rowing.CHROM))}{filesuf}"
        tab = pd.read_csv(filein1, sep='\t')
        left = floor(float(rowing.BPLEFT))
        right = ceil(float(rowing.BPRIGHT))
        tab = tab[(tab['BPWINDOW'] >= left) & (tab['BPWINDOW'] <= right)]
        tab['WEIGHT'] = tab['COUNT'] / tab['COUNT'].sum()
        tab['CUMSUM'] = np.cumsum(tab['WEIGHT'])

        # central tendency
        shape0 = tab.shape[0]
        mx = 0
        moCM = 0
        moBP = 0
        for j in range(shape0):
            row = tab.iloc[j]
            if row['COUNT'] > mx:
                moCM = row['CMWINDOW']
                moBP = row['BPWINDOW']
                mx = row['COUNT']
        moBP = int(moBP)
        meCM = (tab['WEIGHT'] * tab['CMWINDOW']).sum()
        meBP = (tab['WEIGHT'] * tab['BPWINDOW']).sum()
        meBP = int(meBP)
        mdCM = tab[tab['CUMSUM'] >= 0.5]['CMWINDOW'].tolist()[0]
        mdBP = tab[tab['CUMSUM'] >= 0.5]['BPWINDOW'].tolist()[0]
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
    finacol = ['NAME'] + initcol + ['MODEL', 'ALPHA']
    roitab.sort_values(by='MAXIBD', ascending=False, inplace=True)
    nrow = roitab.shape[0]
    roitab['NAME'] = [f'hit{i}' for i in range(1, nrow+1)]
    roitab['MODEL'] = ['a' for _ in range(1, nrow+1)]
    roitab['ALPHA'] = [0.95 for _ in range(1, nrow+1)]
    roitab = roitab[finacol]
    roitab.to_csv(fileout, sep='\t', index=False, columns=['NAME', 'CHROM', 'MINIBD', 'MAXIBD', 'BPCENTER', 'CMCENTER', 'BPLEFTCENTER', 'BPRIGHTCENTER', 'MODEL', 'ALPHA'])

if __name__ == "__main__":
    main()
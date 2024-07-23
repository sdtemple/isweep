# Computation of IBD counts over windows.
# Seth Temple, GitHub: sdtemple
# Last modified: 2024-07-18

import argparse
import numpy as np
import pandas as pd

def main():
    # Set up argparse
    parser = argparse.ArgumentParser(
        description="Compute IBD counts over windows."
    )
    
    # Define arguments
    parser.add_argument(
        'input_file', 
        type=str, 
        help="Input file containing IBD data."
    )
    
    parser.add_argument(
        'output_file', 
        type=str, 
        help="Output file for saving the IBD counts."
    )
    
    parser.add_argument(
        'map_file', 
        type=str, 
        help="Genetic map file."
    )
    
    parser.add_argument(
        '--by', 
        type=int,
        default=20, 
        help="(default: 20) Window size in kb."
    )
    
    parser.add_argument(
        '--end', 
        type=int,
        default=1000, 
        help="(default: 1000) End position in Mb."
    )
    
    parser.add_argument(
        '--telocut', 
        type=float,
        default=2.0, 
        help="(default: 2.0) Drop this much data at telomeres in cM."
    )
    
    # Parse the arguments
    args = parser.parse_args()

    # Type casting and variable setup
    start = 0
    end = int(float(args.end)) * 1_000_000
    by = int(float(args.by)) * 1_000
    cmcut = float(args.telocut)

    # Formatting
    table = pd.read_csv(args.input_file, sep='\t', compression='gzip', header=None)
    columns = list(table.columns)
    cmcol = columns[7]
    encol = columns[6]
    stcol = columns[5]
    table = table[table[encol] >= start]
    table = table[table[stcol] <= end]
    windows = [i for i in range(start, end + by, by)]

    # Setting up for genetic distance
    with open(args.map_file, 'r') as f:
        currline = f.readline().strip()
        currlist = currline.split('\t')
        currcm = float(currlist[2])
        currbp = float(currlist[3])
        prevcm = currcm
        prevbp = currbp
        currline = f.readline().strip()
        currlist = currline.split('\t')
        currcm = float(currlist[2])
        currbp = float(currlist[3])
        cms = []

        # Getting genetic distance
        for w in windows:
            while w >= currbp:
                if currline == '':
                    bpend = currbp
                    cmend = currcm
                    break
                prevline = currline
                prevcm = currcm
                prevbp = currbp
                prevline = currline
                currline = f.readline().strip()
                if currline != '':
                    currlist = currline.split('\t')
                    currcm = float(currlist[2])
                    currbp = float(currlist[3])
                else:
                    currcm = prevcm
                    currbp = prevbp
                    bpend = currbp
                    cmend = currcm
            if currline == '':
                bpend = currbp
                cmend = currcm
            else:
                b = currbp - w
                a = w - prevbp
                d = a + b
                A = a / d
                B = b / d
                cm = B * prevcm + A * currcm
                cms.append(cm)
    
    # Counting IBD segments
    table = table[table[encol] <= bpend]
    windows = [w for w in windows if w <= bpend]
    counts = [((table[stcol] <= i) & (table[encol] >= i)).sum() for i in windows]
    counter = {'BPWINDOW': windows, 'CMWINDOW': cms, 'COUNT': counts}
    counter = pd.DataFrame(counter)

    # CM cutting
    counter = counter[counter['CMWINDOW'] <= (cmend - cmcut)]
    counter = counter[counter['CMWINDOW'] >= cmcut]

    # Count cutting
    counter = counter[counter['COUNT'] > 0]

    # Saving
    counter.to_csv(args.output_file, sep='\t', index=False)

if __name__ == "__main__":
    main()

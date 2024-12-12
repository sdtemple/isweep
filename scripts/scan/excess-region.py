# Combine windows in a region of excess IBD.
# Seth Temple, GitHub: sdtemple
# Last modified: 2024-12-12

import argparse
import pandas as pd

def main():
    # Set up argparse
    parser = argparse.ArgumentParser(
        description="Combine windows in a region of excess IBD."
    )
    
    # Define arguments
    parser.add_argument(
        '--input_file', 
        type=str,
        required=True, 
        help="Input file containing IBD data."
    )
    
    parser.add_argument(
        '--output_file', 
        type=str,
        required=True, 
        help="Output file for saving combined windows of excess IBD."
    )
    
    parser.add_argument(
        '--max_cM_gap', 
        type=float,
        default=0.5, 
        help="(default: 0.5) Allowed cM gap in contiguous stretch of excess IBD."
    )

    parser.add_argument(
        '--statistic', 
        type=str,
        default='COUNT', 
        help="(default: COUNT) Test statistic that exceeds threshold."
    )
    
    # Parse the arguments
    args = parser.parse_args()

    # File handling and computation
    cm = args.max_cM_gap
    g = open(args.output_file, 'w')
    statistic = args.statistic
    minstat = 'MIN' + statistic
    maxstat = 'MAX' + statistic
    table = pd.read_csv(args.input_file,sep='\t')
    nrows = table.shape[0]
    g.write('CHROM\tBPLEFT\tBPRIGHT\tCMLEFT\tCMRIGHT\t' + minstat + '\t' + maxstat + '\n')
    try:
        row = table.iloc[0]
        prev = row
        towrite = prev
        maxct = float(prev[statistic])
        minct = float(prev[statistic])
        for i in range(1,nrows):
            row = table.iloc[i]
            if int(row['CHROM']) == int(prev['CHROM']):
                if float(row['CMWINDOW']) >= (float(prev['CMWINDOW']) + cm):
                    # gap in chromosome
                    g.write(str(int(towrite['CHROM']))); g.write('\t')
                    g.write(str(towrite['BPWINDOW'])); g.write('\t')
                    g.write(str(prev['BPWINDOW'])); g.write('\t')
                    g.write(str(towrite['CMWINDOW'])); g.write('\t')
                    g.write(str(prev['CMWINDOW'])); g.write('\t')
                    g.write(str(minct)); g.write('\t')
                    g.write(str(maxct)); g.write('\n')
                    towrite = row
                    maxct = float(row[statistic])
                    minct = float(row[statistic])
            else:
                # chromosome changed
                g.write(str(int(towrite['CHROM']))); g.write('\t')
                g.write(str(towrite['BPWINDOW'])); g.write('\t')
                g.write(str(prev['BPWINDOW'])); g.write('\t')
                g.write(str(towrite['CMWINDOW'])); g.write('\t')
                g.write(str(prev['CMWINDOW'])); g.write('\t')
                g.write(str(minct)); g.write('\t')
                g.write(str(maxct)); g.write('\n')
                towrite = row
                maxct = float(row[statistic])
                minct = float(row[statistic])
            prev = row
            nextct = float(prev[statistic])
            if nextct > maxct:
                maxct = nextct
            if nextct < minct:
                minct = nextct
        g.write(str(int(towrite['CHROM']))); g.write('\t')
        g.write(str(towrite['BPWINDOW'])); g.write('\t')
        g.write(str(prev['BPWINDOW'])); g.write('\t')
        g.write(str(towrite['CMWINDOW'])); g.write('\t')
        g.write(str(prev['CMWINDOW'])); g.write('\t')
        g.write(str(minct)); g.write('\t')
        g.write(str(maxct)); g.write('\n')
    except IndexError:
        # this happens if there is no excess IBD
        # make a blank file
        pass
    g.close()

if __name__ == "__main__":
    main()

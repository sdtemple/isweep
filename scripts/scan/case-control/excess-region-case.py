# Combine windows in a region of excess IBD.
# Seth Temple, GitHub: sdtemple
# Last modified: 2024-07-22

import argparse

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
    
    # Parse the arguments
    args = parser.parse_args()

    # File handling and computation
    cm = args.max_cM_gap
    g = open(args.output_file, 'w')
    g.write('CHROM\tBPLEFT\tBPRIGHT\tCMLEFT\tCMRIGHT\tMINZDIFFZ\tMAXDIFFZ\n')
    try:
        with open(args.input_file, 'r') as f:
            f.readline()
            line = f.readline()
            row = line.strip().split('\t')
            prev = row
            towrite = prev
            maxct = int(float(prev[7]))
            minct = int(float(prev[7]))
            for line in f:
                row = line.strip().split('\t')
                if int(row[0]) == int(prev[0]):
                    if float(row[1]) >= (float(prev[1]) + cm):
                        # gap in chromosome
                        g.write(towrite[3]); g.write('\t')
                        g.write(towrite[0]); g.write('\t')
                        g.write(prev[0]); g.write('\t')
                        g.write(towrite[1]); g.write('\t')
                        g.write(prev[1]); g.write('\t')
                        g.write(str(minct)); g.write('\t')
                        g.write(str(maxct)); g.write('\n')
                        towrite = row
                        maxct = int(float(row[7]))
                        minct = int(float(row[7]))
                else:
                    # chromosome changed
                    g.write(towrite[3]); g.write('\t')
                    g.write(towrite[0]); g.write('\t')
                    g.write(prev[0]); g.write('\t')
                    g.write(towrite[1]); g.write('\t')
                    g.write(prev[1]); g.write('\t')
                    g.write(str(minct)); g.write('\t')
                    g.write(str(maxct)); g.write('\n')
                    towrite = row
                    maxct = int(float(row[2]))
                    minct = int(float(row[2]))
                prev = row
                nextct = int(float(prev[2]))
                if nextct > maxct:
                    maxct = nextct
                if nextct < minct:
                    minct = nextct
        g.write(towrite[3]); g.write('\t')
        g.write(towrite[0]); g.write('\t')
        g.write(prev[0]); g.write('\t')
        g.write(towrite[1]); g.write('\t')
        g.write(prev[1]); g.write('\t')
        g.write(str(minct)); g.write('\t')
        g.write(str(maxct)); g.write('\n')
    except IndexError:
        # this happens if there is no excess IBD
        # make a blank file
        pass
    g.close()

if __name__ == "__main__":
    main()

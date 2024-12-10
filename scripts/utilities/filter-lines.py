import gzip
import sys
import argparse

parser = argparse.ArgumentParser(description="Apply one filter to tabular data.")

parser.add_argument('--input_file', 
                    type=str, 
                    help='Input tabular data file')
parser.add_argument('--output_file', 
                    type=str, 
                    help='Output tabular data file')
parser.add_argument('--column_index',
                    type=int,
                    default=8,
                    help='Integer [1-indexed] of filtering column'
                    )
parser.add_argument('--lower_bound',
                    type=float,
                    default=0.,
                    help='a in [a,b)'
                    )
parser.add_argument('--upper_bound',
                    type=float,
                    default=2.,
                    help='b in [a,b)'
                    )
parser.add_argument('--complement',
                    type=int,
                    default=1,
                    help='If not 0, then filter is (-inf,a) and [b,inf]'
                    )
parser.add_argument('--comment',
                    type=str,
                    default='I',
                    help='Character denoting lines to skip'
                    )
parser.add_argument('--gzipped',
                    type=int,
                    default=1,
                    help='If 0, the file is not gzipped'
                    )
parser.add_argument('--splitter',
                    type=str,
                    default='\t',
                    help='Character to split lines on'
                    )

args = parser.parse_args()
idx = args.column_index - 1
a = args.lower_bound
b = args.upper_bound
c = args.complement
comment = args.comment
splitter = args.splitter

if args.gzipped:

    # make binary version
    splitter = bytes(splitter, 'utf-8')

    # initalize output
    g = gzip.open(args.output_file,'wb')

    # read input line by line
    with gzip.open(args.input_file,'rb') as f:

        for line in f:
            first_letter = chr(line[0])
            if first_letter == comment:
                g.write(line)

            # complement applied
            elif c:
                full_row = line.strip().split(splitter)
                val = float(full_row[idx])
                if (val < a) or (val >= b):
                    g.write(splitter.join(full_row))
                    g.write(b'\n')
                else:
                    pass 
            
            # complement not applied
            else:
                full_row = line.strip().split(splitter)
                val = float(full_row[idx])
                if (val >= a) and (val < b):
                    g.write(splitter.join(full_row))
                    g.write(b'\n')
                else:
                    pass
            
else:

    # initialize output
    g = open(args.output_file,'w')

    # read input line by line
    with open(args.input_file,'r') as f:
        for line in f:
            if line[0] == comment:
                pass

            # complement applied
            elif c:
                full_row = line.strip().split(splitter)
                val = float(full_row[idx])
                if (val < a) or (val >= b):
                    g.write(splitter.join(full_row))
                    g.write('\n')
                else:
                    pass

            # complement not applied
            else:
                full_row = line.strip().split(splitter)
                val = float(full_row[idx])
                if (val >= a) and (val < b):
                    g.write(splitter.join(full_row))
                    g.write('\n')
                else:
                    pass

g.close()
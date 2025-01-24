import argparse
import numpy as np
import pandas as pd
import multiprocessing as mp
from copy import deepcopy
import functools

def process_chunk(table, map_file_bp, start, end, ind1, ind2, casedict):
    columns = list(table.columns)
    encol = columns[end]
    stcol = columns[start]
    ind1col = columns[ind1]
    ind2col = columns[ind2]
    table[ind1col] = table[ind1col].astype(str)
    table[ind2col] = table[ind2col].astype(str)

    # Map case status to individuals
    table['case1'] = table[ind1col].map(casedict)
    table['case2'] = table[ind2col].map(casedict)
    table['match'] = table['case1'] == table['case2']
    table['casemult'] = table['case1'] * table['case2']
    table['case'] = table.apply(lambda x: x['casemult'] if x['match'] else 2, axis=1)

    # Counting IBD segments
    counts = [((table[stcol] <= i) & (table[encol] >= i)).sum() for i in map_file_bp]
    counts0 = [((table[stcol] <= i) & (table[encol] >= i) & (table['case'] == 0)).sum() for i in map_file_bp]
    counts1 = [((table[stcol] <= i) & (table[encol] >= i) & (table['case'] == 1)).sum() for i in map_file_bp]
    return np.array(counts), np.array(counts0), np.array(counts1)

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description='Compute IBD counts over windows in case and control groups.')
    parser.add_argument('--input_ibd_file', type=str, required=True, help='Input file path')
    parser.add_argument('--input_map_file', type=str, required=True, help='Map file path')
    parser.add_argument('--input_case_file', type=str, required=True, help='File with case status for individuals')
    parser.add_argument('--output_file', type=str, required=True, help='Output file path')
    parser.add_argument('--start', type=int, default=5, help='Start bp column index')
    parser.add_argument('--end', type=int, default=6, help='End bp column index')
    parser.add_argument('--ind1', type=int, default=0, help='Individual 1 column index')
    parser.add_argument('--ind2', type=int, default=2, help='Individual 2 column index')
    parser.add_argument('--chunksize', type=int, default=10000000, help='The parameter in pandas read_csv')

    # Parse the arguments
    args = parser.parse_args()
    input_file = args.input_ibd_file
    map_file = args.input_map_file
    output_file = args.output_file
    start = args.start
    end = args.end
    ind1 = args.ind1
    ind2 = args.ind2
    chunksize = args.chunksize
    casefile = args.input_case_file

    # Formatting the map file
    map_file = pd.read_csv(map_file, sep='\t', header=None)
    map_file.columns = ['chrom', 'rsid', 'cm', 'bp']

    # Initialize the counter
    counts = np.zeros(map_file.shape[0], dtype=int)
    counter = {
        'BPWINDOW': map_file['bp'].to_list(),
        'CMWINDOW': map_file['cm'].to_list(),
        'COUNT': counts,
        'COUNT0': deepcopy(counts),
        'COUNT1': deepcopy(counts),
    }

    # Process case file
    casedict = dict()
    with open(casefile) as f:
        for line in f:
            ind, status = line.strip().split('\t')
            casedict[ind] = int(float(status))

    # Set up a pool of worker processes
    pool = mp.Pool(mp.cpu_count())

    try:
        process_func = functools.partial(process_chunk, map_file_bp=map_file['bp'], start=start, end=end, ind1=ind1, ind2=ind2, casedict=casedict)
        results = []

        for chunk in pd.read_csv(input_file, sep='\t', chunksize=chunksize):
            result = pool.apply_async(process_func, args=(chunk,))
            results.append(result)

        pool.close()
        pool.join()

        for result in results:
            counts_chunk, counts0_chunk, counts1_chunk = result.get()
            counter['COUNT'] += counts_chunk
            counter['COUNT0'] += counts0_chunk
            counter['COUNT1'] += counts1_chunk

    finally:
        pool.terminate()

    counter_df = pd.DataFrame(counter)
    # Count cutting
    counter_df = counter_df[counter_df['COUNT'] > 0]

    # Saving
    counter_df.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    main()
# Computing IBD counts over windows

# importing
import argparse
import numpy as np
import pandas as pd
import multiprocessing as mp

# Function to process each chunk
def process_chunk(chunk, start, end, map_file_bp):
    columns = list(chunk.columns)
    encol = columns[end]
    stcol = columns[start]
    counts = [((chunk[stcol] <= i) & (chunk[encol] >= i)).sum() for i in map_file_bp]
    return counts

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description='Compute IBD counts over windows.')
    parser.add_argument('--input_ibd_file', 
                        type=str,
                        required=True, 
                        help='Input file path')
    parser.add_argument('--input_map_file', 
                        type=str,
                        required=True, 
                        help='Map file path')
    parser.add_argument('--output_file', 
                        type=str,
                        required=True, 
                        help='Output file path')
    parser.add_argument('--start', 
                        type=int,
                        default=5, 
                        help='Start bp column index')
    parser.add_argument('--end', 
                        type=int,
                        default=6, 
                        help='End bp column index')
    parser.add_argument('--chunksize',
                        type=int,
                        default=10000000,
                        help='The parameter in pandas read_csv'
                        )

    # Parse the arguments
    args = parser.parse_args()
    input_file = args.input_ibd_file
    map_file = args.input_map_file
    output_file = args.output_file
    start = args.start
    end = args.end
    chunksize = args.chunksize

    # formatting the map file
    map_file = pd.read_csv(map_file, sep='\t', header=None)
    map_file.columns = ['chrom', 'rsid', 'cm', 'bp']

    # initialize the counter
    counts = np.zeros(map_file.shape[0], dtype=int)
    counter = {
        'BPWINDOW': map_file['bp'].to_list(),
        'CMWINDOW': map_file['cm'].to_list(),
        'COUNT': counts
    }

    # Setup multiprocessing pool
    pool = mp.Pool(mp.cpu_count())

    results = []
    try:
        # counting IBD segments in chunks
        for chunk in pd.read_csv(input_file, sep='\t', chunksize=chunksize):
            result = pool.apply_async(process_chunk, args=(chunk, start, end, map_file['bp']))
            results.append(result)

        pool.close()  # Close the pool to new tasks
        pool.join()   # Wait for all worker processes to exit

        # Combine results from all workers
        for result in results:
            counts_chunk = np.array(result.get())
            counter['COUNT'] += counts_chunk

    finally:
        # Make sure to close the pool in case of an error
        pool.terminate()

    counter = pd.DataFrame(counter)

    # count cutting
    counter = counter[counter['COUNT'] > 0]

    # saving
    counter.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    main()
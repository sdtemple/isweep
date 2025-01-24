import gzip
import argparse
import multiprocessing as mp
from functools import partial

def process_lines(lines, idx, a, b, c, splitter, comment):
    results = []
    for line in lines:
        if line.startswith(comment):
            results.append(line)
        else:
            full_row = line.strip().split(splitter)
            val = float(full_row[idx])
            if (c and ((val < a) or (val >= b))) or (not c and (a <= val < b)):
                results.append(line)
    return results

def write_output(output_file, is_gzipped, all_results):
    mode = 'wb' if is_gzipped else 'w'
    with (gzip.open(output_file, mode) if is_gzipped else open(output_file, mode)) as f_out:
        for result in all_results:
            for line in result:
                f_out.write(line.encode('utf-8') if is_gzipped else line)

def read_in_chunks(file_obj, chunk_size=1024*1024):
    while True:
        data = file_obj.readlines(chunk_size)
        if not data:
            break
        yield data

def filter_file(input_file, output_file, idx, a, b, c, comment, splitter, gzipped, chunksize):
    # Define the function to process chunks of lines
    process_func = partial(process_lines, idx=idx, a=a, b=b, c=c, splitter=splitter, comment=comment)
    
    if gzipped:
        with gzip.open(input_file, 'rt', encoding='utf-8') as f_in:
            chunks = list(read_in_chunks(f_in, chunksize))
    else:
        with open(input_file, 'r', encoding='utf-8') as f_in:
            chunks = list(read_in_chunks(f_in, chunksize))
    
    # Set up a pool of worker processes
    with mp.Pool(mp.cpu_count()) as pool:
        all_results = pool.map(process_func, chunks)

    # Write results to output file
    write_output(output_file, gzipped, all_results)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Apply one filter to tabular data.")
    parser.add_argument('--input_file', type=str, help='Input tabular data file')
    parser.add_argument('--output_file', type=str, help='Output tabular data file')
    parser.add_argument('--column_index', type=int, default=8, help='Integer [1-indexed] of filtering column')
    parser.add_argument('--lower_bound', type=float, default=0., help='Lower bound of filter range [a, b)')
    parser.add_argument('--upper_bound', type=float, default=2., help='Upper bound of filter range [a, b)')
    parser.add_argument('--complement', type=int, default=1, help='Complement filter range: (-inf, a) and [b, inf]')
    parser.add_argument('--comment', type=str, default='I', help='Character denoting lines to skip')
    parser.add_argument('--gzipped', type=int, default=1, help='Is the file gzipped? 0 = No, 1 = Yes')
    parser.add_argument('--splitter', type=str, default='\t', help='Character to split lines on')
    parser.add_argument('--chunksize', type=int, default=10000000, help='Number of lines per process')

    args = parser.parse_args()
    idx = args.column_index - 1
    a = args.lower_bound
    b = args.upper_bound
    c = args.complement
    comment = args.comment
    splitter = args.splitter

    filter_file(args.input_file, args.output_file, idx, a, b, c, comment, splitter, args.gzipped, args.chunksize)
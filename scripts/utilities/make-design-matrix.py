import os
import pandas as pd
import argparse

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process phenotype and outlier data to make design matrix.")
    parser.add_argument("--output_design_matrix_prefix", type=str, help="Prefix for the output design matrix files.")
    parser.add_argument("--input_phenotype_file", type=str, help="Path to the input phenotype file.")
    parser.add_argument("--input_outlier_prefix", type=str, default='outlier', help="Prefix for the input outlier files.")
    parser.add_argument("--input_outlier_suffix", type=str, default='.tsv', help="Suffix for the input outlier files.")
    parser.add_argument("--first_index", type=int, default=1, help="The starting index of outlier files.")

    # Parse arguments
    args = parser.parse_args()

    input_phenotype_file = args.input_phenotype_file
    input_prefix = args.input_outlier_prefix
    input_suffix = args.input_outlier_suffix
    first_idx = args.first_index
    output_design_matrix_prefix = args.output_design_matrix_prefix

    itr = first_idx
    while os.path.exists(input_prefix + str(itr) + input_suffix):
        itr += 1

    dict_phenotype = {}
    with open(input_phenotype_file, 'r') as f:
        for line in f:
            line = line.strip().split()
            dict_phenotype[str(line[0]) + '_1'] = line[1]
            dict_phenotype[str(line[0]) + '_2'] = line[1]

    with open(output_design_matrix_prefix + '.tsv', 'w') as out_phenotype:
        indicator_columns = ['OUTLIER' + str(j) for j in range(first_idx, itr)]
        out_phenotype.write('\t'.join(['INDIVIDUAL', 'HAPLOTYPE', 'PHENO'] + indicator_columns + ['OUTLIER_ANY']) + '\n')
        outliers = []
        for i in range(first_idx, itr):
            with open(input_prefix + str(i) + input_suffix, 'r') as f:
                for line in f:
                    line = line.strip().split()
                    out_phenotype.write('\t'.join([str(line[0])[:-2], str(line[0]), line[1]]))
                    for j in range(i - 1):
                        out_phenotype.write('\t0')
                    out_phenotype.write('\t1')
                    for j in range(i + 1, itr):
                        out_phenotype.write('\t0')
                    out_phenotype.write('\t1\n')
                    outliers.append(str(line[0]))
        less_phenotypes = {k: v for k, v in dict_phenotype.items() if v not in outliers}
        for k, v in less_phenotypes.items():
            out_phenotype.write('\t'.join([str(k)[:-2], str(k), str(v)]))
            for j in range(itr):
                out_phenotype.write('\t0')
            out_phenotype.write('\n')

    table = pd.read_csv(output_design_matrix_prefix + '.tsv', sep='\t')
    table.sort_values(by=['INDIVIDUAL', 'HAPLOTYPE'], inplace=True)
    table.to_csv(output_design_matrix_prefix + '.sorted.tsv', sep='\t', index=False)

if __name__ == "__main__":
    main()
# Write per-outlier IBD subset tables from a headerless hap-ibd file.
# Seth Temple, GitHub: sdtemple
# Added outlier IBD subset workflow support on 2026-04-20

import argparse
import os
import pandas as pd


def read_outlier_groups(input_outlier_prefix, input_outlier_suffix, first_index):
    """Read sequential outlier group files and return (index, haplotype-set) pairs."""
    groups = []
    idx = first_index
    while os.path.exists(input_outlier_prefix + str(idx) + input_outlier_suffix):
        path = input_outlier_prefix + str(idx) + input_outlier_suffix
        haplotypes = []
        with open(path, 'r') as f:
            for line in f:
                row = line.strip().split('\t')
                if len(row) == 0 or row[0] == '':
                    continue
                haplotypes.append(str(row[0]))
        groups.append((idx, set(haplotypes)))
        idx += 1
    return groups


def main():
    parser = argparse.ArgumentParser(
        description='Write one .ibd.gz subset per outlier file using haplotype membership.'
    )
    parser.add_argument(
        '--input_ibd_file',
        type=str,
        required=True,
        help='Path to headerless hap-ibd file (e.g., narrowing.filt.case.ibd.gz).',
    )
    parser.add_argument(
        '--input_outlier_prefix',
        type=str,
        required=True,
        help='Prefix for outlier files (e.g., cohort/hit/outlier).',
    )
    parser.add_argument(
        '--input_outlier_suffix',
        type=str,
        default='.phenotype.tsv',
        help='Suffix for outlier files (default: .phenotype.tsv).',
    )
    parser.add_argument(
        '--first_index',
        type=int,
        default=1,
        help='First outlier index to scan (default: 1).',
    )
    parser.add_argument(
        '--input_phenotype_file',
        type=str,
        required=True,
        help='Two-column phenotype file: sample ID and phenotype.',
    )
    parser.add_argument(
        '--output_suffix_file',
        type=str,
        default='.ibd.gz',
        help='Suffix for output subset files (default: .ibd.gz).',
    )
    parser.add_argument(
        '--output_ibd_prefix',
        type=str,
        required=True,
        help='Prefix for output subset files (e.g., cohort/hit/outlier).',
    )
    args = parser.parse_args()

    groups = read_outlier_groups(
        args.input_outlier_prefix,
        args.input_outlier_suffix,
        args.first_index,
    )

    ibd = pd.read_csv(
        args.input_ibd_file,
        sep='\t',
        header=None,
        compression='infer',
        usecols=[0, 1, 2, 3, 4, 5, 6, 7],
        names=['IND1', 'HAP1', 'IND2', 'HAP2', 'CHROM', 'LEFT', 'RIGHT', 'CM'],
    )

    phen = pd.read_csv(
        args.input_phenotype_file,
        sep='\t',
        header=None,
        usecols=[0, 1],
        names=['IND', 'PHENO'],
    )
    pheno_map = dict(zip(phen['IND'].astype(str), phen['PHENO']))

    ibd['PHENO1'] = ibd['IND1'].astype(str).map(pheno_map)
    ibd['PHENO2'] = ibd['IND2'].astype(str).map(pheno_map)

    ibd['PAIR1'] = ibd['IND1'].astype(str) + '_' + ibd['HAP1'].astype(str)
    ibd['PAIR2'] = ibd['IND2'].astype(str) + '_' + ibd['HAP2'].astype(str)

    for idx, group_haplotypes in groups:
        subset = ibd[
            ibd['PAIR1'].isin(group_haplotypes)
            & ibd['PAIR2'].isin(group_haplotypes)
        ][['IND1', 'HAP1', 'IND2', 'HAP2', 'CHROM', 'LEFT', 'RIGHT', 'CM', 'PHENO1', 'PHENO2']]

        output_file = (
            args.output_ibd_prefix
            + str(idx)
            + args.output_suffix_file
        )
        subset.to_csv(output_file, sep='\t', index=False, compression='gzip')


if __name__ == '__main__':
    main()

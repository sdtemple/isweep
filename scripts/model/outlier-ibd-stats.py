# Summarize IBD segment endpoints for outlier haplotype groups.
# Seth Temple, GitHub: sdtemple
# Added endpoint summary workflow support on 2026-04-20

import argparse
import os
import pandas as pd


def read_outlier_groups(input_outlier_prefix, input_outlier_suffix, first_index):
    """Read sequential outlier group files and return group labels with haplotypes."""
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
        groups.append((f'outlier{idx}', set(haplotypes)))
        idx += 1
    return groups


def summarize_group_segments(ibd_table, group_haplotypes):
    """Compute endpoint summary statistics for segments where both haps are in group."""
    default_chrom = pd.NA
    if ibd_table.shape[0] > 0:
        default_chrom = ibd_table['CHROM'].mode().iloc[0]

    if len(group_haplotypes) == 0:
        return {
            'chromosome': default_chrom,
            'n_haplotypes': 0,
            'n_segments': 0,
            'left_median': pd.NA,
            'left_q75': pd.NA,
            'left_q90': pd.NA,
            'left_max': pd.NA,
            'right_min': pd.NA,
            'right_q10': pd.NA,
            'right_q25': pd.NA,
            'right_median': pd.NA,
        }

    sub = ibd_table[
        ibd_table['HAP1'].isin(group_haplotypes)
        & ibd_table['HAP2'].isin(group_haplotypes)
    ]

    if sub.shape[0] == 0:
        return {
            'chromosome': default_chrom,
            'n_haplotypes': len(group_haplotypes),
            'n_segments': 0,
            'left_median': pd.NA,
            'left_q75': pd.NA,
            'left_q90': pd.NA,
            'left_max': pd.NA,
            'right_min': pd.NA,
            'right_q10': pd.NA,
            'right_q25': pd.NA,
            'right_median': pd.NA,
        }

    left = sub['LEFT_BP']
    right = sub['RIGHT_BP']

    return {
        'chromosome': sub['CHROM'].mode().iloc[0],
        'n_haplotypes': len(group_haplotypes),
        'n_segments': int(sub.shape[0]),
        'left_median': float(left.median()),
        'left_q75': float(left.quantile(0.75)),
        'left_q90': float(left.quantile(0.90)),
        'left_max': float(left.max()),
        'right_min': float(right.min()),
        'right_q10': float(right.quantile(0.10)),
        'right_q25': float(right.quantile(0.25)),
        'right_median': float(right.median()),
    }


def main():
    parser = argparse.ArgumentParser(
        description='Summarize left/right endpoint statistics for outlier IBD segments.'
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
        '--output_file',
        type=str,
        required=True,
        help='Output TSV file with one row per outlier group.',
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
        usecols=[0, 1, 2, 3, 4, 5, 6],
        names=['ID1', 'HAPIDX1', 'ID2', 'HAPIDX2', 'CHROM', 'LEFT_BP', 'RIGHT_BP'],
    )

    # This workflow expects one ROI per run, so IBD rows must come from one chromosome.
    chrom_values = pd.unique(ibd['CHROM'].dropna())
    if len(chrom_values) > 1:
        raise ValueError(
            'Expected IBD table with one chromosome, but found multiple: '
            + ', '.join([str(c) for c in sorted(chrom_values)])
        )

    ibd['HAP1'] = ibd['ID1'].astype(str) + '_' + ibd['HAPIDX1'].astype(str)
    ibd['HAP2'] = ibd['ID2'].astype(str) + '_' + ibd['HAPIDX2'].astype(str)

    rows = []
    for group_name, group_haplotypes in groups:
        stats = summarize_group_segments(ibd, group_haplotypes)
        stats['outlier_group'] = group_name
        rows.append(stats)

    out = pd.DataFrame(rows)
    if out.shape[0] == 0:
        out = pd.DataFrame(columns=[
            'outlier_group',
            'chromosome',
            'n_haplotypes',
            'n_segments',
            'left_median',
            'left_q75',
            'left_q90',
            'left_max',
            'right_min',
            'right_q10',
            'right_q25',
            'right_median',
        ])
    else:
        out = out[
            [
                'outlier_group',
                'chromosome',
                'n_haplotypes',
                'n_segments',
                'left_median',
                'left_q75',
                'left_q90',
                'left_max',
                'right_min',
                'right_q10',
                'right_q25',
                'right_median',
            ]
        ]

    out.to_csv(args.output_file, sep='\t', index=False)


if __name__ == '__main__':
    main()

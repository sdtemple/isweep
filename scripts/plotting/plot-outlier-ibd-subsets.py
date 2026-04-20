# Plot sampled IBD segments from per-outlier IBD subset files.
# Seth Temple, GitHub: sdtemple
# Added outlier IBD subset plotting support on 2026-04-20

import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter


def list_outlier_ibd_files(input_prefix, input_suffix, first_index):
    """List sequential outlier IBD files and return (index, file_path) pairs."""
    files = []
    idx = first_index
    while os.path.exists(input_prefix + str(idx) + input_suffix):
        files.append((idx, input_prefix + str(idx) + input_suffix))
        idx += 1
    return files


def read_locus_center(locus_file):
    """Read CENTER value from locus.case.txt style key-value file."""
    with open(locus_file, 'r') as f:
        for line in f:
            row = line.strip().split('\t')
            if len(row) < 2:
                continue
            if row[0] == 'CENTER':
                return float(row[1])
    raise ValueError(f'CENTER not found in locus file: {locus_file}')


def main():
    parser = argparse.ArgumentParser(
        description='Plot sampled segments from outlierN.ibd.gz files.'
    )
    parser.add_argument(
        '--input_ibd_prefix',
        type=str,
        required=True,
        help='Prefix for input IBD subset files (e.g., cohort/hit/outlier).',
    )
    parser.add_argument(
        '--input_ibd_suffix',
        type=str,
        default='.ibd.gz',
        help='Suffix for input IBD subset files (default: .ibd.gz).',
    )
    parser.add_argument(
        '--first_index',
        type=int,
        default=1,
        help='First outlier index to scan (default: 1).',
    )
    parser.add_argument(
        '--output_plot_prefix',
        type=str,
        required=True,
        help='Prefix for output plot files (e.g., cohort/hit/outlier).',
    )
    parser.add_argument(
        '--output_plot_suffix',
        type=str,
        default='.ibd.png',
        help='Suffix for output plot files (default: .ibd.png).',
    )
    parser.add_argument(
        '--sample_size',
        type=int,
        default=50,
        help='(default: 50) Number of segments to sample without replacement.',
    )
    parser.add_argument(
        '--random_seed',
        type=int,
        default=20260420,
        help='(default: 20260420) Random seed used for sampling.',
    )
    parser.add_argument(
        '--width',
        type=float,
        default=10.0,
        help='(default: 10.0) Figure width in inches.',
    )
    parser.add_argument(
        '--height',
        type=float,
        default=5.0,
        help='(default: 5.0) Figure height in inches.',
    )
    parser.add_argument(
        '--padding_fraction',
        type=float,
        default=0.15,
        help='(default: 0.10) Fractional x-axis padding applied to both sides.',
    )
    parser.add_argument(
        '--segment_color',
        type=str,
        default='black',
        help='(default: black) Color of segment lines.',
    )
    parser.add_argument(
        '--segment_linewidth',
        type=float,
        default=1.0,
        help='(default: 1.0) Line width for segment lines.',
    )
    parser.add_argument(
        '--label_fontsize',
        type=float,
        default=12.0,
        help='(default: 12.0) Font size for x/y labels and title.',
    )
    parser.add_argument(
        '--input_locus_file',
        type=str,
        required=True,
        help='Path to locus.case.txt file containing CENTER.',
    )
    args = parser.parse_args()

    center = read_locus_center(args.input_locus_file)

    files = list_outlier_ibd_files(
        args.input_ibd_prefix,
        args.input_ibd_suffix,
        args.first_index,
    )

    for idx, file_path in files:
        table = pd.read_csv(file_path, sep='\t', compression='infer')
        if table.shape[0] == 0:
            continue

        pheno_counts_text = ''
        if all(col in table.columns for col in ['IND1', 'IND2', 'PHENO1', 'PHENO2']):
            individuals = pd.concat(
                [
                    table[['IND1', 'PHENO1']].rename({'IND1': 'IND', 'PHENO1': 'PHENO'}, axis=1),
                    table[['IND2', 'PHENO2']].rename({'IND2': 'IND', 'PHENO2': 'PHENO'}, axis=1),
                ],
                ignore_index=True,
            )
            individuals = individuals.drop_duplicates(subset=['IND'])
            counts = individuals['PHENO'].value_counts(dropna=False)
            lines = ['Phenotype']
            for k, v in counts.items():
                lines.append(f'{k}: {int(v)}')
            pheno_counts_text = '\n'.join(lines)

        chrom_values = pd.unique(table['CHROM'].dropna())
        if len(chrom_values) != 1:
            raise ValueError(
                f'Expected one chromosome in {file_path}, found: '
                + ', '.join([str(c) for c in sorted(chrom_values)])
            )
        chrom = str(chrom_values[0])

        n = min(args.sample_size, table.shape[0])
        total_segments = int(table.shape[0])
        sampled = table.sample(n=n, replace=False, random_state=args.random_seed).copy()
        sampled.reset_index(drop=True, inplace=True)

        left_min = float(sampled['LEFT'].min())
        right_max = float(sampled['RIGHT'].max())
        span = right_max - left_min
        pad = args.padding_fraction * span
        if span == 0:
            pad = max(1.0, abs(left_min) * args.padding_fraction)

        y = [i + 1 for i in range(sampled.shape[0])]

        plt.figure(figsize=(args.width, args.height))
        plt.hlines(
            y=y,
            xmin=sampled['LEFT'],
            xmax=sampled['RIGHT'],
            colors=args.segment_color,
            linewidth=args.segment_linewidth,
        )
        plt.axvline(center, color='tab:red', linestyle='--', linewidth=1.25)
        plt.title(
            f'Outlier {idx} ({n} sampled out of {total_segments} segments)',
            fontsize=args.label_fontsize,
        )
        plt.xlabel(f'Chromosome {chrom}', fontsize=args.label_fontsize)
        plt.ylabel('Haplotype pairs', fontsize=args.label_fontsize)
        plt.yticks([])
        plt.xlim(left_min - pad, right_max + pad)

        ax = plt.gca()
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((0, 0))
        ax.xaxis.set_major_formatter(formatter)
        ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
        ax.grid(True, axis='x', color='lightgray', linestyle='-', linewidth=0.6, alpha=0.75)
        if pheno_counts_text:
            ax.text(
                0.01,
                0.98,
                pheno_counts_text,
                transform=ax.transAxes,
                va='top',
                ha='left',
                fontsize=max(8.0, args.label_fontsize - 2.0),
                bbox=dict(facecolor='white', edgecolor='lightgray', alpha=1),
            )

        output_file = args.output_plot_prefix + str(idx) + args.output_plot_suffix
        plt.tight_layout()
        plt.savefig(output_file, dpi=300)
        plt.close()


if __name__ == '__main__':
    main()

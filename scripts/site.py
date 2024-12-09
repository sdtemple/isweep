# Refine sweep site

# imports
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Set up the argument parser
parser = argparse.ArgumentParser(description='Refine sweep site from SNP data.')

parser.add_argument('--input_snp_file', 
                    type=str,
                    required=True, 
                    help='Path to the input SNP file (table 1)')
parser.add_argument('--output_folder', 
                    type=str,
                    required=True, 
                    help='Folder for output')
parser.add_argument('--window_index', 
                    type=int,
                    default=0, 
                    help='Column index for distance unit')
parser.add_argument('--freq_index', 
                    type=int,
                    default=1, 
                    help='Column index for frequency')
parser.add_argument('--window_size', 
                    type=float, 
                    help='Size of sliding window')
parser.add_argument('--window_step', 
                    type=float, 
                    help='Step of sliding window')
parser.add_argument('--quantile_sum', 
                    type=int, 
                    help='Sum 1... quantiles')
parser.add_argument('--max_spacing', 
                    type=float, 
                    help='Max spacing filter')

# Parse the arguments
args = parser.parse_args()
snpin = args.input_snp_file
folderout = args.output_folder
winidx = args.window_index
freqidx = args.freq_index
winsize = args.window_size
winstep = args.window_step
qrng = args.quantile_sum
maxspace = args.max_spacing

# functions
def q_score(table,
            freqidx=1,
            winidx=0,
            by=0.05,
            size=0.2,
            qrange=10,
            jump=0.01
           ):
    '''Compute quantile score over sliding window

    Parameters
    ----------
    table : pandas DataFrame
        SNPs with position, allele frequency
    freqidx : int
        Frequency column index
    winidx : int
        Distance column index
    by : int
        Step for sliding window
    size : numeric
        Size for sliding window
    qrange : int
        Sum quantiles 1...
    jump : int
        Filter if spacing exceeds this in a window

    Returns
    -------
    pandas DataFrame
        Table for quantile scored windows
    '''
    # column names
    headernames = list(table.columns)
    pos = headernames[winidx]  # position name
    aaf = headernames[freqidx]  # frequency name
    # spacings
    table.sort_values(by=pos, ignore_index=True, inplace=True)
    spacings = []
    x = table[pos]
    for i in range(table.shape[0]):
        try:
            spacings.append(x[i + 1] - x[i])
        except KeyError:  # for last item
            spacings.append(0)
    table['SPACING'] = spacings
    table['JUMP'] = table['SPACING'] >= jump
    # sliding windows
    winmin = table[pos].min()
    winmax = table[pos].max() - size
    winrange = np.arange(winmin, winmax, by)
    windowed = dict()
    ctr = 0
    for w in winrange:
        winleft = w
        winright = w + size
        winmid = (winleft + winright) / 2
        # boundaries
        subtable = table[(table[pos] >= winleft) & (table[pos] < winright)]
        # has a SNP
        sz = subtable.shape[0]
        if sz >= 1:
            # has no jump
            if subtable['JUMP'].sum() <= 0:
                sm = 0
                for i in range(qrange):
                    per = (i + 1) * 0.01
                    val = np.quantile(subtable[aaf], per)
                    sm += val
                sp = subtable['SPACING'].max()
                windowed[ctr] = (winmid,
                                 sm,
                                 sz,
                                 sp
                                 )
                ctr += 1
    # data frame
    windowed = pd.DataFrame(windowed)
    windowed = windowed.T
    windowed.columns = ['POS', 'QSCORE', 'SIZE', 'MAXSPACING']
    windowed.sort_values('QSCORE', ascending=False, inplace=True)
    return windowed

# compute quantile score over sliding window
snptab = pd.read_csv(snpin, sep='\t')
qtab = q_score(snptab,
               freqidx,
               winidx,
               winstep,
               winsize,
               qrng,
               maxspace
               )
qtab.to_csv(f"{folderout}/first.qs.tsv.gz",
            sep='\t',
            index=False
            )

# best window
best = list(qtab.iloc[0])
bestbp = best[0]

# plotting
qtab.sort_values(by='POS', ascending=True, inplace=True)
plt.plot(qtab['POS'],
         qtab['QSCORE'],
         color='k'
         )
plt.ylabel('Score')
plt.xlabel('Position')
plt.savefig(f"{folderout}/first.score.png")
plt.clf()

headernames = list(snptab.columns)
wincol = headernames[winidx]  # position name
freqcol = headernames[freqidx]  # frequency name
plt.scatter(snptab[wincol],
            snptab[freqcol],
            color='k',
            s=5
            )
plt.ylim(-0.1, 1.1)
plt.ylabel('SNP frequency')
plt.xlabel('Position')
plt.savefig(f"{folderout}/first.snp.png", dpi=300)
plt.clf()

# record positions
with open(f"{folderout}/first.pos.txt", 'w') as f:
    f.write('bp\t')
    f.write(str(int(bestbp)))
    f.write('\n')
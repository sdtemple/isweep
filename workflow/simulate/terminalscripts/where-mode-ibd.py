import sys
import gzip
input_file, = sys.argv[1:]
f = gzip.open(input_file, 'rt')
f.readline()
windows = []
counts = []
for line in f:
    bpwindow, cmwindow, count = line.strip().split('\t')
    windows.append(int(float(bpwindow)))
    counts.append(int(float(count)))
max_val = max(counts)
idx = counts.index(max_val)
# caution: this will find the first window
# w/ the mode ibd segment count
print(windows[idx])
f.close()

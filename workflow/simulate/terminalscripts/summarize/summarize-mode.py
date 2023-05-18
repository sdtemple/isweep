
import os
import sys
import gzip
import logging

# inputting
dirin, filein, fileout = sys.argv[1:]

g = open(dirin + '/' + fileout, 'w')

modepos=[]
modepos2=[]
for subdir in os.listdir(dirin):
    input_file = dirin + '/' + subdir + '/' + filein
    try:
        f = gzip.open(input_file, 'rt')
        f.readline()
        windows = []
        counts = []
        for line in f:
            bpwindow, cmwindow, count = line.strip().split('\t')
            # windows.append(int(float(bpwindow)))
            windows.append(float(cmwindow))
            # use comments for cM based analysis
            counts.append(int(float(count)))
        max_val = max(counts)
        idx = counts.index(max_val)
        counts.reverse()
        idx2 = counts.index(max_val)
        idx2 = len(counts) - idx2
        if idx != idx2:
            logging.warning('There are at least 2 modes.')
        modepos.append(windows[idx])
        modepos2.append(windows[idx2])
        # caution: this will find the first window
        # when processing from left and from right
        # w/ the mode ibd segment count
        f.close()
    except NotADirectoryError:
        pass

# print a summary
avg=sum(modepos)/len(modepos)
avg2=sum(modepos2)/len(modepos2)
print('Average of left mode: '+str(avg))
print('Average of right mode: '+str(avg2))
minmode=min(modepos)
maxmode=max(modepos2)
print('Minimum of left mode :'+str(minmode))
print('Maximum of right mode :'+str(maxmode))

# writing mode file over replicate study
g.write('MODELEFT\tMODERIGHT\n')
num=len(modepos)
for n in range(num):
    m=modepos[n]
    m2=modepos2[n]
    g.write(str(m)); g.write('\t')
    g.write(str(m2)); g.write('\n')
g.close()

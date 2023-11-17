import sys
folder, fileout, percent_rule, plurality_rule, sample_size = sys.argv[1:]
percent_rule = float(percent_rule)
plurality_rule =float(plurality_rule)
sample_size = float(sample_size)
ctr= 0 
numnodes = []
try:
    while True:
        ctr += 1
        numnode = 0
        with open(folder + '/outlier' + str(ctr) + '.txt', 'r') as f:
            for line in f:
                numnode += 1
        numnodes.append(numnode)
except:
    pass
g=open(fileout,'w')
largest_percent = max(numnodes) / sample_size
plurality_percent = max(numnodes) / sum(numnodes)
if (largest_percent >= percent_rule) and (plurality_percent >= plurality_rule):
    g.write('1\n')
    print('1')
else:
    g.write('0\n')
    print('0')
g.close()
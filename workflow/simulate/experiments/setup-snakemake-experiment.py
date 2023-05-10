import sys
K,start,micro,s,t,fileout=sys.argv[1:]
K=int(float(K)) # num replicates
start=int(float(start)) # index to begin at for sim name
s=float(s) # selection coefficient
t=int(float(t)) # time of adaptive mutation
with open(fileout,'w') as f:
    f.write('MICROEXP\tSIMNAME\tSELCOEF\tTIMEMUT\n')
    for k in range(K):
        idx=k+start
        f.write(micro); f.write('\t')
        f.write(str(idx)); f.write('\t')
        f.write(str(s)); f.write('\t')
        f.write(str(t)); f.write('\n')

# # write another one for pop struct
# with open(fileout,'w') as f:
#     f.write('MICROEXP\tSIMNAME\tSELCOEF\tTIMEMUT\n')
#     for k in range(K):
#         idx=k+start
#         f.write(micro); f.write('\t')
#         f.write(str(idx)); f.write('\t')
#         f.write(str(s)); f.write('\t')
#         f.write(str(t)); f.write('\n')

# output the outlier haplotypes with their phenotypes
import sys
input_outlier_file,input_phenotype_file,output_file,=sys.argv[1:]
phen=dict()
with open(input_phenotype_file,'r') as f:
    for line in f:
        hap,sta=line.strip().split('\t')
        phen[hap]=int(sta)
print(phen)
outl=dict()
with open(input_outlier_file,'r') as f:
    for line in f:
        print(line.strip())
        hap=line.strip()
        outl[hap]=phen[hap[:-2]]
print(outl)
print(outl.items())
with open(output_file,'w') as g:
    #kv=[k,v for k,v in outl.items()]
    for k,v in outl.items():
        g.write(str(k))
        g.write('\t')
        g.write(str(v))
        g.write('\n')

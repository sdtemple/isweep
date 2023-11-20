"""
Make uniform map file (cM distances)

@author: sdtemple
"""

rho = float(snakemake.config['CHANGE']['SIMULATE']['RHO'])
L = float(snakemake.config['CHANGE']['SIMULATE']['CMLEN'])
out=snakemake.output.mapout
rho = float(rho) * 100
L = int(float(L))
cM = L
RHO = float(rho)
L = float(L)
L = L / RHO
L = int(L) - 1
with open(out, 'w') as f:
    f.write('1\t.\t0.0\t0\n')
    f.write('1\t.\t' + str(cM) + '\t' + str(L) + '\n')

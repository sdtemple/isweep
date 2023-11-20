import sys
out, rho, L = sys.argv[1:]
rho = float(rho)
L = float(L)
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

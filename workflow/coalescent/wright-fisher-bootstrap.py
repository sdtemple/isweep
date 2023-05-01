#!/bin/bash

# running wright-fisherian walks for bootstrap inference
# seth d temple, sdtemple@uw.edu
# may 1, 2023

import sys
from isweep import *

filein,fileout,s,p,Ne,inh,tau,sv=sys.argv[1:]
s=float(s)
p=float(p)
tau=int(float(tau))
sv=float(sv)

tab=pd.read_csv()

g=open(fileout,'w')
g.write('TRUES','TRUEPARAM','PARAM','')

import sys

outfolder, tsvfile, numpcs, outRscript = sys.argv[1:]
numpcs=int(float(numpcs))+1

script='#!/bin/usr/env Rscript\ntable=read.table("'
script+=tsvfile
script+='",header=T,sep="\\t")\n'
for i in range(1,numpcs):
    for j in range(i+1,numpcs):
        script+='png("'
        script+=outfolder
        # main
        script+='/pcaplot.'
        script+=str(i)
        script+='.'
        script+=str(j)
        script+='.png")\n'
        script+='plot(table$EV'
        script+=str(i)
        script+=',table$EV'
        script+=str(j)
        script+=',xlab=\'Eigenvector '
        script+=str(i)
        script+='\',ylab=\'Eigenvector '
        script+=str(j)
        script+='\',cex.lab=1.5,cex.axis=1.5,cex.main=1.5,cex.sub=1.5)\n'
f=open(outRscript,'w')
f.write(script)
f.close()

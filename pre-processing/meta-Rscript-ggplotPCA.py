import sys

outfolder, tsvfile, numpcs, covariate, outRscript = sys.argv[1:]
numpcs=int(float(numpcs))+1

script='#!/bin/usr/env Rscript\nlibrary(ggplot2)\ntable=read.table("'
script+=tsvfile
script+='",header=T,sep="\\t")\n'
script+='table$'+covariate+'=factor(table$'+covariate+')\n'
for i in range(1,numpcs):
    for j in range(i+1,numpcs):
        script+='png("'
        script+=outfolder
        script+='/pcaplot.'
        script+=str(i)
        script+='.'
        script+=str(j)
        script+='.'
        script+=covariate
        script+='.png")\n'
        script+='p=ggplot(table,aes(EV'
        script+=str(i)
        script+=',EV'
        script+=str(j)
        script+=',color='
        script+=covariate
        script+='))+xlab(\'Eigenvector '
        script+=str(i)
        script+='\')+ylab(\'Eigenvector '
        script+=str(j)
        script+='\')+geom_point(shape=1)+theme_classic()\n'
        script+='print(p, vp=grid::viewport(gp=grid::gpar(cex=1.5)))\n'
        script+='dev.off()\n'
        script+='\n'
f=open(outRscript,'w')
f.write(script)
f.close()

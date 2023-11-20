import sys

outfolder, tsvfile, numpcs, covariate, downsample, outRscript = sys.argv[1:]
numpcs=int(float(numpcs))+1

script='#!/bin/usr/env Rscript\nlibrary(ggplot2)\ntable=read.table("'
script+=tsvfile
script+='",header=T,sep="\\t")\n'
script+='table$'+covariate+'=factor(table$'+covariate+')\n'
script+='table=table[table$'+covariate+'!= 0,]\n'
script+='table=table[!is.na(table$'+covariate+'),]\n'
script+='facs=unique(table$'+covariate+')\n'
script+='keep=c()\n'
script+='for(fac in facs){\n\tindices=which(table$'+covariate+'==fac,arr.ind=T)\n\tif(length(indices)<'+downsample+'){\n\t\tkeep=c(keep,indices)\n\t} else{\n\t\tkeep=c(keep,sample(indices,'+downsample+'))\n\t}\n}\n'
script+='table=table[keep,]\n\n'
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

# Manhattan IBD Counts

# setup, imports ----------------------------------------------------------

# importing
library(ggplot2)
library(dplyr)

# inputting
args=commandArgs(trailingOnly=T)
suffix=args[2] # .ibd.windowed.tsv.gz
# align with config file
folder=args[1] # [CHANGE][FOLDERS][STUDY]
chrlow=args[3] # [CHANGE][ISWEEP][CHRLOW]
chrhigh=args[4] # [CHANGE][ISWEEP][CHRLOW]
cutoff1=args[5] # [FIXED][ISWEEP][SCANCUTOFF]
cutoff2=args[6] # [FIXED][ISWEEP][TELOCUTOFF]

# type casting
chrlow=as.integer(chrlow)
chrhigh=as.integer(chrhigh)

# reading in data
ibd = read.table(paste(folder,'/ibdsegs/ibdends/modified/scan/','chr',chrlow,suffix, sep = ''), header = T, stringsAsFactors = F)
ibd$CHROM = chrlow
for(i in (chrlow+1):(chrhigh)){
  print(i)
  incoming = read.table(paste(folder,'/ibdsegs/ibdends/modified/scan/','chr',i,suffix, sep = ''), header = T, stringsAsFactors = F)
  incoming$CHROM = i
  ibd = rbind(ibd, incoming)
}


# maths -------------------------------------------------------------------

# remove some obvious outliers (telo/centromeres, ROIs)
medi=median(ibd$COUNT)
stdv=sd(ibd$COUNT)
a=medi-stdv*cutoff2
b=medi+stdv*cutoff2
low=a
newibd=ibd[ibd$COUNT>=a,]
newibd=newibd[ibd$COUNT<=b,]

# compute browning and browning 2020 cutoff
medi=median(newibd$COUNT)
stdv=sd(newibd$COUNT)
upp=medi+stdv*cutoff1
md=medi

# chromosome maths
ibd$CUMPOS = NA
s = 0
mbp = c()
chr = unique(ibd$CHROM)
nchr = length(chr)
for(i in chrlow:chrhigh){
  mbp[i] = max(ibd[ibd$CHROM == chr[i],]$CMWINDOW)
  ibd[ibd$CHROM == chr[i], "CUMPOS"] = ibd[ibd$CHROM == chr[i], "CMWINDOW"] + s
  s = s + mbp[i]
}
axis = ibd %>% group_by(CHROM) %>% summarize(center = (max(CUMPOS) + min(CUMPOS)) / 2)

# plotting ----------------------------------------------------------------

# colors
gb = rep(c('gray','black'), times = ceil((chrhigh-chrlow)/2))
ibd$CHROM = as.numeric(ibd$CHROM)

# counts
ggplot(ibd, aes(x = CUMPOS, y = COUNT, color = as.factor(CHROM))) +
  geom_point(size = 0.1) +
  scale_x_continuous(label = axis$CHROM, breaks = axis$center) +
  scale_color_manual(values = gb) +
  geom_abline(slope = 0, intercept = upp, color = 'red') +
  geom_abline(slope = 0, intercept = low, color = 'orange') +
  geom_abline(slope = 0, intercept = md, color = 'pink') +
  labs(x = NULL, y = 'IBD Segment Count') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position = 'none')

# saving plots
typ=paste('.manhattan.count.','high',cutoff1,'.low',cutoff2,'.',sep='')
dp='print'
pic='jpeg'
ggsave(paste(folder,'/plots/','chr',chrlow,'-',chrhigh,typ,pic,sep=''),device=pic,dpi=dp)
pic='png'
ggsave(paste(folder,'/plots/','chr',chrlow,'-',chrhigh,typ,pic,sep=''),device=pic,dpi=dp)
pic='tiff'
ggsave(paste(folder,'/plots/','chr',chrlow,'-',chrhigh,typ,pic,sep=''),device=pic,dpi=dp)
pic='eps'
ggsave(paste(folder,'/plots/','chr',chrlow,'-',chrhigh,typ,pic,sep=''),device=pic,dpi=dp)

# log 10 counts
ibd$LOG10COUNT = log10(ibd$COUNT+1)
ggplot(ibd, aes(x = CUMPOS, y = LOG10COUNT, color = as.factor(CHROM))) +
  geom_point(size = 0.1) +
  scale_x_continuous(label = axis$CHROM, breaks = axis$center) +
  scale_color_manual(values = gb) +
  geom_abline(slope = 0, intercept = log10(upp), color = 'red') +
  geom_abline(slope = 0, intercept = log10(low), color = 'orange') +
  geom_abline(slope = 0, intercept = log10(md), color = 'pink') +
  labs(x = NULL, y = 'Log10 IBD Segment Count') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position = 'none')

# saving plots
typ=paste('.manhattan.log10count.','high',cutoff1,'.low',cutoff2,'.',sep='')
dp='print'
pic='jpeg'
ggsave(paste(folder,'/plots/','chr',chrlow,'-',chrhigh,typ,pic,sep=''),device=pic,dpi=dp)
pic='png'
ggsave(paste(folder,'/plots/','chr',chrlow,'-',chrhigh,typ,pic,sep=''),device=pic,dpi=dp)
pic='tiff'
ggsave(paste(folder,'/plots/','chr',chrlow,'-',chrhigh,typ,pic,sep=''),device=pic,dpi=dp)
pic='eps'
ggsave(paste(folder,'/plots/','chr',chrlow,'-',chrhigh,typ,pic,sep=''),device=pic,dpi=dp)

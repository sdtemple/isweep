# Implement the OU process -based GW significance threshold for IBD rate
# Seth D Temple
# May 24, 2024

start_time <- Sys.time()

# Processing Inputs -------------------------------------------------------

args=commandArgs(trailingOnly=T)
fileout = args[1] # name of output file
theta = args[2] # ou process correlation parameter
pval = args[3] # significance threshold
stepsize = args[4] # step size in Morgans (Delta in S&Y[2007] )
chrnum = args[5] # number of chromosomes in Morgans
chrlen = args[6] # length of each chromosome in Morgans
J = args[7] # to determine the simulation-based threshold
M = args[8] # for calculating the type 1 error 

the.lines = paste('theta: ',theta,'\n',sep='')
the.lines = paste(the.lines,'significance-threshold: ',pval,'\n',sep='')
the.lines = paste(the.lines,'step-size: ',stepsize,'\n',sep='')
the.lines = paste(the.lines,'number-of-chromosomes: ',chrnum,'\n',sep='')
the.lines = paste(the.lines,'length-of-chromosomes: ',chrlen,'\n',sep='')
the.lines = paste(the.lines,'amount-of-sim-training: ',J,'\n',sep='')
the.lines = paste(the.lines,'amount-to-compute-type1: ',M,'\n',sep='')

theta = as.numeric(theta)
pval = as.numeric(pval)
siglevel = pval
stepsize = as.numeric(stepsize)
chrnum = as.integer(chrnum)
chrlen = as.integer(chrlen)
J = as.integer(J)
M = as.integer(M)


# Local functions ---------------------------------------------------------


# from siegmund and yakir (2007) the statistics of gene mapping
Nu = function(y){
  y = y / 2
  (1 / y) * (pnorm(y) - 0.5) / (y * pnorm(y) + dnorm(y))
}

# from siegmund and yakir (2007) the statistics of gene mapping
OU.approx = function(z,beta,Delta,length,chr,center,test="one-sided"){
  d = switch(test,"two-sided"=2,"one-sided"=1)
  p = 1 - exp(-d*chr*(1-pnorm(z))-d*beta*length*z*dnorm(z)*Nu(z*sqrt(2*beta*Delta)))
  return(p-center)
}

# from siegmund and yakir (2007) the statistics of gene mapping
OU.approx.cont = function(z,beta,length,chr,center,test="one-sided"){
  d = switch(test,"two-sided"=2,"one-sided"=1)
  p = 1 - exp(-d*chr*(1-pnorm(z))-d*beta*length*z*dnorm(z))
  return(p-center)
}


# Training the simulation-based approach ----------------------------------


K=floor(chrlen/stepsize)*chrnum

{
  maxs = c()
  for(j in 1:J){
    if((j%%100)==0){print(j)}
    x=rnorm(1)
    xs = c(x)
    for(i in 2:K){
      newx=rnorm(1,mean=x*exp(-stepsize*theta),sd=sqrt(2-2*exp(-stepsize*theta)))
      xs = c(xs,newx)
      x=newx
    }
    zs=(xs-mean(xs))/sd(xs)
    maxs = c(maxs,max(zs))
  }
}

print("simulation base finished")

# Calculating Type 1 error ------------------------------------------------

ctr.cont = 0
ctr.disc = 0
ctr.simu = 0
ctr.bonf = 0

lin = uniroot(OU.approx.cont,c(1,20),theta,chrnum*chrlen,chrnum,siglevel)$root
lin2 = uniroot(OU.approx,c(1,20),theta,stepsize,chrnum*chrlen,chrnum,siglevel)$root
lin3 = quantile(maxs, 1 - pval)
lin4 = qnorm(1-pval/K)

print("the thresholds are:")
print(c(lin,lin2,as.numeric(lin3),lin4))

{
  for(m in 1:M){
    if((m%%100)==0){print(m)}
    x=rnorm(1)
    xs = c(x)
    for(i in 2:K){
      newx=rnorm(1,mean=x*exp(-stepsize*theta),sd=sqrt(2-2*exp(-stepsize*theta)))
      xs = c(xs,newx)
      x=newx
    }
    zs=(xs-mean(xs))/sd(xs)
    if(length(zs[zs>=lin2])>0){ctr.disc=ctr.disc+1}
    if(length(zs[zs>=lin])>0){ctr.cont=ctr.cont+1}
    if(length(zs[zs>=lin3])>0){ctr.simu=ctr.simu+1}
    if(length(zs[zs>=lin4])>0){ctr.bonf=ctr.bonf+1}
  }
}

# Writing output ----------------------------------------------------------

fileConn<-file(fileout)

the.lines = paste(the.lines,'quantile-discrete-approx: ',lin2,'\n',sep='')
the.lines = paste(the.lines,'quantile-continuous-approx: ',lin,'\n',sep='')
the.lines = paste(the.lines,'quantile-simulation-based: ',lin3,'\n',sep='')
the.lines = paste(the.lines,'quantile-bonferroni: ',lin4,'\n',sep='')
the.lines = paste(the.lines,'type1-discrete-approx: ',ctr.disc/M,'\n',sep='')
the.lines = paste(the.lines,'type1-continuous-approx: ',ctr.cont/M,'\n',sep='')
the.lines = paste(the.lines,'type1-simulation-based: ',ctr.simu/M,'\n',sep='')
the.lines = paste(the.lines,'type1-bonferroni: ',ctr.bonf/M,'\n',sep='')

end_time <- Sys.time()

elapsed_time <- end_time - start_time

the.lines = paste(the.lines,'system-time: ',elapsed_time,sep='')

writeLines(the.lines,fileConn)

close(fileConn)

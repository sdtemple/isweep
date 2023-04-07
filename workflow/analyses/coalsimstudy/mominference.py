import sys
import os
from isweep import *

# first arg: # replicates
# second arg: # bootstraps
# third arg: # diploids
# fourth arg: param fixed (s, p, or Ne)
# fifth arg: param fixed (s, p, or Ne)
# sixth arg: example bins file
# seventh arg: output folder
# eigth arg: output file prefix

folder=sys.argv[7]
if not os.path.exists(folder):
    os.mkdir(folder)

# fixed parameter settings
# change selection coefficients if interested
lst=[0.02,0.03,0.04] # selection coefficients
# lst=[0.02,0.03] # selection coefficients

# input parameter settings
p=float(sys.argv[4]) # allele freq
Ne=read_Ne(sys.argv[5]) # demo history
bn=read_bins(sys.argv[6])

# input sim settings
nreplicates=int(float(sys.argv[1]))
nbootstraps=int(float(sys.argv[2]))
nsamples=int(float(sys.argv[3]))

# fixed sim settings
ploidy=2
msamples=ploidy*nsamples
long_ibd=3.0
N=msamples*(msamples-1)/2-msamples
outfile=folder+'/'+sys.argv[8]

# write to file concurrently

f=open(outfile+'.good.tsv','w')
f.write('PARAM\tTRUE\tINITEST\tCORREST\tCORRLOW\tCORRUPP\n')

# goodness of fit to length distribution
bn.append(np.inf)
ab=bn
sinfs=[[] for k in range(len(lst))]
for i in range(nreplicates):
    for k in range(len(lst)):
        s=lst[k]
        out = simulate_ibd_isweep(
            nsamples,
            s,
            p,
            Ne,
            long_ibd,
            long_ibd,
            ploidy=ploidy
        )
        ibd=big_format_distribution(out[0][2],out[0][4])
        ibd=bin_ibd_segments(ibd,ab)
        se = minimize_scalar(
            chi2_isweep,
            args=(p,Ne,N,ibd,ab),
            bounds=(0,0.5),
            method='bounded'
        ).x
        sboot=[]
        for j in range(nbootstraps):
            print(lst[k],i,j) # comment out to stop stdout
            out = simulate_ibd_isweep(
                nsamples,
                s,
                p,
                Ne,
                long_ibd,
                long_ibd,
                ploidy=ploidy
            )
            ibd=big_format_distribution(out[0][2],out[0][4])
            ibd=np.array(bin_ibd_segments(ibd,ab))
            sj = minimize_scalar(
                chi2_isweep,
                args=(p,Ne,N,ibd,ab),
                bounds=(0,0.5),
                method='bounded'
            ).x
            sboot.append(sj)
        sint=bootstrap_standard_bc(se,sboot)
        sinf=['GOOD',s,se,sint[1],sint[0],sint[2]]
        sinfs[k].append(sinf)
        # saving
        row=sinf
        for l in range(len(row)):
            f.write(str(row[l]))
            if l != (len(row)-1):
                f.write('\t')
            else:
                f.write('\n')

f.close()

# write to file concurrently

f=open(outfile+'.lbls.tsv','w')
f.write('PARAM\tTRUE\tINITEST\tCORREST\tCORRLOW\tCORRUPP\n')

# goodness of fit to length distribution (labeled)
# bn.append(np.inf)
# ab=bn
sinfs=[[] for k in range(len(lst))]
for i in range(nreplicates):
    for k in range(len(lst)):
        s=lst[k]
        out = simulate_ibd_isweep(
            nsamples,
            s,
            p,
            Ne,
            long_ibd,
            long_ibd,
            ploidy=ploidy
        )
        ibd0=big_format_distribution(out[2][2],out[2][4])
        ibd0=bin_ibd_segments(ibd0,ab)
        ibd1=big_format_distribution(out[1][2],out[1][4])
        ibd1=bin_ibd_segments(ibd1,ab)
        se = minimize_scalar(
            chi2_labeled_isweep,
            args=(p,Ne,N,ibd1,ibd0,ab),
            bounds=(0,0.5),
            method='bounded'
        ).x
        sboot=[]
        for j in range(nbootstraps):
            print(lst[k],i,j) # comment out to stop stdout
            out = simulate_ibd_isweep(
                nsamples,
                s,
                p,
                Ne,
                long_ibd,
                long_ibd,
                ploidy=ploidy
            )
            ibd0=big_format_distribution(out[2][2],out[2][4])
            ibd0=bin_ibd_segments(ibd0,ab)
            ibd1=big_format_distribution(out[1][2],out[1][4])
            ibd1=bin_ibd_segments(ibd1,ab)
            sj = minimize_scalar(
                chi2_labeled_isweep,
                args=(p,Ne,N,ibd1,ibd0,ab),
                bounds=(0,0.5),
                method='bounded'
            ).x
            sboot.append(sj)
        sint=bootstrap_standard_bc(se,sboot)
        sinf=['LBLS',s,se,sint[1],sint[0],sint[2]]
        sinfs[k].append(sinf)
        # saving
        row=sinf
        for l in range(len(row)):
            f.write(str(row[l]))
            if l != (len(row)-1):
                f.write('\t')
            else:
                f.write('\n')

f.close()

# write to file concurrently

f=open(outfile+'.mome.tsv','w')
f.write('PARAM\tTRUE\tINITEST\tCORREST\tCORRLOW\tCORRUPP\n')

# method of moments to ibd count
ab=[long_ibd,np.inf]
sinfs=[[] for k in range(len(lst))]
for i in range(nreplicates):
    for k in range(len(lst)):
        s=lst[k]
        out = simulate_ibd_isweep(
            nsamples,
            s,
            p,
            Ne,
            long_ibd,
            long_ibd,
            ploidy=ploidy
        )
        ibd=out[0][0]
        se = minimize_scalar(
            chi2_isweep,
            args=(p,Ne,N,(ibd,),ab),
            bounds=(0,0.5),
            method='bounded'
        ).x
        sboot=[]
        for j in range(nbootstraps):
            print(lst[k],i,j) # comment out to stop stdout
            ibd = simulate_ibd_isweep(
                nsamples,
                se,
                p,
                Ne,
                long_ibd,
                long_ibd,
                ploidy=ploidy
            )
            ibd=out[0][0]
            sj = minimize_scalar(
                chi2_isweep,
                args=(p,Ne,N,(ibd,),ab),
                bounds=(0,0.5),
                method='bounded'
            ).x
            sboot.append(sj)
        sint=bootstrap_standard_bc(se,sboot)
        sinf=['MOME',s,se,sint[1],sint[0],sint[2]]
        sinfs[k].append(sinf)
        # saving
        row=sinf
        for l in range(len(row)):
            f.write(str(row[l]))
            if l != (len(row)-1):
                f.write('\t')
            else:
                f.write('\n')

f.close()

# # write to file at the end
#
# # goodness of fit to length distribution
# bn.append(np.inf)
# ab=bn
# sinfs=[[] for k in range(len(lst))]
# for i in range(nreplicates):
#     for k in range(len(lst)):
#         s=lst[k]
#         out = simulate_ibd_isweep(
#             nsamples,
#             s,
#             p,
#             Ne,
#             long_ibd,
#             long_ibd,
#             ploidy=ploidy
#         )
#         ibd=big_format_distribution(out[0][2],out[0][4])
#         ibd=bin_ibd_segments(ibd,ab)
#         se = minimize_scalar(
#             chi2_isweep,
#             args=(p,Ne,N,ibd,ab),
#             bounds=(0,0.5),
#             method='bounded'
#         ).x
#         sboot=[]
#         for j in range(nbootstraps):
#             print(lst[k],i,j)
#             out = simulate_ibd_isweep(
#                 nsamples,
#                 s,
#                 p,
#                 Ne,
#                 long_ibd,
#                 long_ibd,
#                 ploidy=ploidy
#             )
#             ibd=big_format_distribution(out[0][2],out[0][4])
#             ibd=np.array(bin_ibd_segments(ibd,ab))
#             sj = minimize_scalar(
#                 chi2_isweep,
#                 args=(p,Ne,N,ibd,ab),
#                 bounds=(0,0.5),
#                 method='bounded'
#             ).x
#             sboot.append(sj)
#         sint=bootstrap_standard_bc(se,sboot)
#         sinf=['GOOD',s,se,sint[1],sint[0],sint[2]]
#         sinfs[k].append(sinf)
#
# with open(outfile+'.good.tsv','w') as f:
#     f.write('PARAM\tTRUE\tINITEST\tCORREST\tCORRLOW\tCORRUPP\n')
#     for i in range(len(sinfs)):
#         sinf=sinfs[i]
#         for j in range(len(sinf)):
#             row=sinf[j]
#             for k in range(len(row)):
#                 f.write(str(row[k]))
#                 if k != (len(row)-1):
#                     f.write('\t')
#                 else:
#                     f.write('\n')
#
# # goodness of fit to length distribution (labeled)
# bn.append(np.inf)
# ab=bn
# sinfs=[[] for k in range(len(lst))]
# for i in range(nreplicates):
#     for k in range(len(lst)):
#         s=lst[k]
#         out = simulate_ibd_isweep(
#             nsamples,
#             s,
#             p,
#             Ne,
#             long_ibd,
#             long_ibd,
#             ploidy=ploidy
#         )
#         ibd0=big_format_distribution(out[2][2],out[2][4])
#         ibd0=bin_ibd_segments(ibd0,ab)
#         ibd1=big_format_distribution(out[1][2],out[1][4])
#         ibd1=bin_ibd_segments(ibd1,ab)
#         se = minimize_scalar(
#             chi2_labeled_isweep,
#             args=(p,Ne,N,ibd1,ibd0,ab),
#             bounds=(0,0.5),
#             method='bounded'
#         ).x
#         sboot=[]
#         for j in range(nbootstraps):
#             print(lst[k],i,j)
#             out = simulate_ibd_isweep(
#                 nsamples,
#                 s,
#                 p,
#                 Ne,
#                 long_ibd,
#                 long_ibd,
#                 ploidy=ploidy
#             )
#             ibd0=big_format_distribution(out[2][2],out[2][4])
#             ibd0=bin_ibd_segments(ibd0,ab)
#             ibd1=big_format_distribution(out[1][2],out[1][4])
#             ibd1=bin_ibd_segments(ibd1,ab)
#             sj = minimize_scalar(
#                 chi2_labeled_isweep,
#                 args=(p,Ne,N,ibd1,ibd0,ab),
#                 bounds=(0,0.5),
#                 method='bounded'
#             ).x
#             sboot.append(sj)
#         sint=bootstrap_standard_bc(se,sboot)
#         sinf=['LBLS',s,se,sint[1],sint[0],sint[2]]
#         sinfs[k].append(sinf)
#
# with open(outfile+'.lbls.tsv','w') as f:
#     f.write('PARAM\tTRUE\tINITEST\tCORREST\tCORRLOW\tCORRUPP\n')
#     for i in range(len(sinfs)):
#         sinf=sinfs[i]
#         for j in range(len(sinf)):
#             row=sinf[j]
#             for k in range(len(row)):
#                 f.write(str(row[k]))
#                 if k != (len(row)-1):
#                     f.write('\t')
#                 else:
#                     f.write('\n')
#
#
# # method of moments to ibd count
# ab=[long_ibd,np.inf]
# sinfs=[[] for k in range(len(lst))]
# for i in range(nreplicates):
#     for k in range(len(lst)):
#         s=lst[k]
#         out = simulate_ibd_isweep(
#             nsamples,
#             s,
#             p,
#             Ne,
#             long_ibd,
#             long_ibd,
#             ploidy=ploidy
#         )
#         ibd=out[0][0]
#         se = minimize_scalar(
#             chi2_isweep,
#             args=(p,Ne,N,(ibd,),ab),
#             bounds=(0,0.5),
#             method='bounded'
#         ).x
#         sboot=[]
#         for j in range(nbootstraps):
#             print(lst[k],i,j)
#             ibd = simulate_ibd_isweep(
#                 nsamples,
#                 se,
#                 p,
#                 Ne,
#                 long_ibd,
#                 long_ibd,
#                 ploidy=ploidy
#             )
#             ibd=out[0][0]
#             sj = minimize_scalar(
#                 chi2_isweep,
#                 args=(p,Ne,N,(ibd,),ab),
#                 bounds=(0,0.5),
#                 method='bounded'
#             ).x
#             sboot.append(sj)
#         sint=bootstrap_standard_bc(se,sboot)
#         sinf=['MOME',s,se,sint[1],sint[0],sint[2]]
#         sinfs[k].append(sinf)
#
# with open(outfile+'.mome.tsv','w') as f:
#     f.write('PARAM\tTRUE\tINITEST\tCORREST\tCORRLOW\tCORRUPP\n')
#     for i in range(len(sinfs)):
#         sinf=sinfs[i]
#         for j in range(len(sinf)):
#             row=sinf[j]
#             for k in range(len(row)):
#                 f.write(str(row[k]))
#                 if k != (len(row)-1):
#                     f.write('\t')
#                 else:
#                     f.write('\n')

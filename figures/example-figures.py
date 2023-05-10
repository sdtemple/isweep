from isweep import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

s=0.03
p=0.70
Ne=read_Ne('bottleneck-1000G.ne')
K=10

# standard with inset
fig, ax = plt.subplots()
axins = ax.inset_axes([0.5, 0.5, 0.47, 0.47])
out=walk_variant_backward(s,p,Ne,one_step_model='a')
plt.plot(range(len(out[0])),out[0],linewidth=1)
axins.plot(range(len(out[0])),np.log10(out[0]),linewidth=1,color="tab:blue")
plt.ylim(-0.1,1.1)
axins.set_ylim(-5,np.log10(1.1))
plt.xlim(-10,500)
axins.set_xlim(-10,500)
plt.ylabel('Allele frequency')
axins.set_ylabel('Log10')
plt.xlabel('Generation')
axins.set_xlabel('Generation')
plt.title('Example selective sweep')
for k in range(K):
    out=walk_variant_backward(s,p,Ne,random_walk=True,one_step_model='a')
    axins.plot(range(len(out[0])),np.log10(out[0]),linewidth=1,color='tab:blue')
    plt.plot(range(len(out[0])),out[0],linewidth=1,color='tab:blue')
fig.savefig('example.sweep.inset.png')

# different allele frequency
fig, ax = plt.subplots()
axins = ax.inset_axes([0.5, 0.5, 0.47, 0.47])
out=walk_variant_backward(s,p,Ne,one_step_model='a')
plt.plot(range(len(out[0])),out[0],linewidth=1)
axins.plot(range(len(out[0])),np.log10(out[0]),linewidth=1,color="tab:blue")
out=walk_variant_backward(s,p+0.1,Ne,one_step_model='a')
plt.plot(range(len(out[0])),out[0],linewidth=1,color='tab:orange')
axins.plot(range(len(out[0])),np.log10(out[0]),linewidth=1,color="tab:orange")
plt.ylim(-0.1,1.1)
axins.set_ylim(-5,np.log10(1.1))
plt.xlim(-10,500)
axins.set_xlim(-10,500)
plt.ylabel('Allele frequency')
axins.set_ylabel('Log10')
plt.xlabel('Generation')
axins.set_xlabel('Generation')
plt.title('Example: different allele frequency')
for k in range(K):
    out=walk_variant_backward(s,p,Ne,random_walk=True,one_step_model='a')
    axins.plot(range(len(out[0])),np.log10(out[0]),linewidth=1,color='tab:blue')
    plt.plot(range(len(out[0])),out[0],linewidth=1,color='tab:blue')
    out=walk_variant_backward(s,p+0.1,Ne,one_step_model='a',random_walk=True)
    plt.plot(range(len(out[0])),out[0],linewidth=1,color='tab:orange')
    axins.plot(range(len(out[0])),np.log10(out[0]),linewidth=1,color="tab:orange")
fig.savefig('example.sweep.pplus.png')

# ongoing sweep or not
fig, ax = plt.subplots()
axins = ax.inset_axes([0.5, 0.5, 0.47, 0.47])
out=walk_variant_backward(s,p,Ne,one_step_model='a')
plt.plot(range(len(out[0])),out[0],linewidth=1)
axins.plot(range(len(out[0])),np.log10(out[0]),linewidth=1,color="tab:blue")
out=walk_variant_backward(s,p,Ne,one_step_model='a',tau0=20)
plt.plot(range(len(out[0])),out[0],linewidth=1,color='tab:orange')
axins.plot(range(len(out[0])),np.log10(out[0]),linewidth=1,color="tab:orange")
plt.ylim(-0.1,1.1)
axins.set_ylim(-5,np.log10(1.1))
plt.xlim(-10,500)
axins.set_xlim(-10,500)
plt.ylabel('Allele frequency')
axins.set_ylabel('Log10')
plt.xlabel('Generation')
axins.set_xlabel('Generation')
plt.title('Example: different allele frequency')
for k in range(K):
    out=walk_variant_backward(s,p,Ne,random_walk=True,one_step_model='a')
    axins.plot(range(len(out[0])),np.log10(out[0]),linewidth=1,color='tab:blue')
    plt.plot(range(len(out[0])),out[0],linewidth=1,color='tab:blue')
    out=walk_variant_backward(s,p,Ne,one_step_model='a',random_walk=True,tau0=20)
    plt.plot(range(len(out[0])),out[0],linewidth=1,color='tab:orange')
    axins.plot(range(len(out[0])),np.log10(out[0]),linewidth=1,color="tab:orange")
fig.savefig('example.sweep.tauplus.png')

# different genetic model
p2=0.3
fig, ax = plt.subplots()
axins = ax.inset_axes([0.5, 0.5, 0.47, 0.47])
# small p
out=walk_variant_backward(s,p2,Ne,one_step_model='a')
plt.plot(range(len(out[0])),out[0],linewidth=1)
axins.plot(range(len(out[0])),np.log10(out[0]),linewidth=1,color="tab:blue")
out=walk_variant_backward(s,p2,Ne,one_step_model='d')
plt.plot(range(len(out[0])),out[0],linewidth=1,color='tab:orange')
axins.plot(range(len(out[0])),np.log10(out[0]),linewidth=1,color="tab:orange")
# large p
out=walk_variant_backward(s,p,Ne,one_step_model='a')
plt.plot(range(len(out[0])),out[0],linewidth=1,color='tab:green')
axins.plot(range(len(out[0])),np.log10(out[0]),linewidth=1,color="tab:green")
out=walk_variant_backward(s,p,Ne,one_step_model='d')
plt.plot(range(len(out[0])),out[0],linewidth=1,color='tab:red')
axins.plot(range(len(out[0])),np.log10(out[0]),linewidth=1,color="tab:red")
plt.ylim(-0.1,1.1)
axins.set_ylim(-5,np.log10(1.1))
plt.xlim(-10,500)
axins.set_xlim(-10,500)
plt.ylabel('Allele frequency')
axins.set_ylabel('Log10')
plt.xlabel('Generation')
axins.set_xlabel('Generation')
plt.title('Example: different genetic model')
for k in range(K):
    # small p
    out=walk_variant_backward(s,p2,Ne,random_walk=True,one_step_model='a')
    axins.plot(range(len(out[0])),np.log10(out[0]),linewidth=1,color='tab:blue')
    plt.plot(range(len(out[0])),out[0],linewidth=1,color='tab:blue')
    out=walk_variant_backward(s,p2,Ne,one_step_model='d',random_walk=True)
    plt.plot(range(len(out[0])),out[0],linewidth=1,color='tab:orange')
    axins.plot(range(len(out[0])),np.log10(out[0]),linewidth=1,color="tab:orange")
    # large p
    out=walk_variant_backward(s,p,Ne,random_walk=True,one_step_model='a')
    axins.plot(range(len(out[0])),np.log10(out[0]),linewidth=1,color='tab:green')
    plt.plot(range(len(out[0])),out[0],linewidth=1,color='tab:green')
    out=walk_variant_backward(s,p,Ne,one_step_model='d',random_walk=True)
    plt.plot(range(len(out[0])),out[0],linewidth=1,color='tab:red')
    axins.plot(range(len(out[0])),np.log10(out[0]),linewidth=1,color="tab:red")
fig.savefig('example.sweep.genetic.png')

# distribution of segment counts by frequency
tab=pd.read_csv('sim.tract.distr.tsv',sep='\t')
sns.kdeplot(tab['0.25'])
sns.kdeplot(tab['0.5'])
sns.kdeplot(tab['0.75'])
plt.xlim(0,5000)
plt.title('Distributions of segment counts')
plt.legend(labels=['0.25','0.50','0.75'],title='Frequency')
plt.savefig('kde.frequency.png')

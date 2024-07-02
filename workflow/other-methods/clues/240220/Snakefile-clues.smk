# clues v1 analysis pipeline

# setup macro folder
import os
macro=str(config['CHANGE']['FOLDERS']['MACRO'])
if not os.path.exists(macro):
	os.mkdir(macro)

# load in experiments, set up micro, simname folders
import pandas as pd
micro=str(config['CHANGE']["FOLDERS"]["MICRO"])
sims = pd.read_csv(micro, sep='\t', header=0)
J = sims.shape[0]
for j in range(J):
	row = sims.loc[j,]
	if not os.path.exists(macro+'/'+str(row.MICROEXP)):
		os.mkdir(macro+'/'+str(row.MICROEXP))
	if not os.path.exists(macro+'/'+str(row.MICROEXP)+'/'+str(row.SIMNAME)):
		os.mkdir(macro+'/'+str(row.MICROEXP)+'/'+str(row.SIMNAME))
sims['FOLDER'] = [macro + '/' + sims.loc[j].MICROEXP + '/' + str(sims.loc[j].SIMNAME) for j in range(J)]
sims['FILE'] = sims['FOLDER'] + '/results.hap.tsv'
sims['EXISTS'] = sims['FILE'].apply(os.path.isfile)
sims=sims[sims['EXISTS']]
sims = sims.set_index("SIMNAME", drop=False)

rule all:
	input:
		[f"{sim.FOLDER}/clues.mcmc.selcoef" for sim in sims.itertuples()],
		[f"{sim.FOLDER}/relate.mcmc.chr1.vcf.gz" for sim in sims.itertuples()],

##### relate ---> clues #####

# format coal file beforehand

# downsample, subset to region near causal mutation
rule relate_vcf:
    input:
        vcfin='{macro}/{micro}/{seed}/second.chr1.vcf.gz',
        subsample=macro+'/'+str(config['CHANGE']['RELATE']['SAMPLE']),
    output:
        vcfout='{macro}/{micro}/{seed}/relate.mcmc.chr1.vcf',
    params:
        loc=str(config['FIXED']['SIMULATE']['LOC']),
        pm=str(config['CHANGE']['RELATE']['BUFFER']),
    shell:
        """
        left=$(python -c "out = {params.loc} - {params.pm} ; print(out)")
        right=$(python -c "out = {params.loc} + {params.pm} ; print(out)")
        tabix -fp vcf {input.vcfin}
        bcftools view {input.vcfin} \
            -r 1:${{left}}-${{right}} \
            -S {input.subsample}  \
            -Ov -o {output.vcfout}
        """

# set up input files from vcf
rule relate_input:
    input:
        vcfin='{macro}/{micro}/{seed}/relate.mcmc.chr1.vcf',
    output:
        hapsout='{macro}/{micro}/{seed}/relate.mcmc.haps',
        sampout='{macro}/{micro}/{seed}/relate.mcmc.sample',
    params:
        relate=str(config['CHANGE']['FOLDERS']['SOFTWARE'])+'/'+str(config['CHANGE']['PROGRAMS']['RELATE']),
        prefix='{macro}/{micro}/{seed}/relate.mcmc.chr1',
    shell:
        """
        {params.relate}/bin/RelateFileFormats \
            --mode ConvertFromVcf \
            --haps {output.hapsout} \
            --sample {output.sampout} \
            --input {params.prefix}
        """

# run relate for first time
rule relate_run:
    input:
        hapsout='{macro}/{micro}/{seed}/relate.mcmc.haps',
        sampout='{macro}/{micro}/{seed}/relate.mcmc.sample',
        genmap='{macro}/uniform.map',
        coa=macro+'/'+str(config['CHANGE']['RELATE']['COAL']),
    output:
        anc='{macro}/{micro}/{seed}/relate.mcmc.anc',
        mut='{macro}/{micro}/{seed}/relate.mcmc.mut',
    params:
        relate=str(config['CHANGE']['FOLDERS']['SOFTWARE'])+'/'+str(config['CHANGE']['PROGRAMS']['RELATE']),
        mu=str(config['FIXED']['SIMULATE']['MU']),
        mid='{micro}.{seed}',
    shell:
        """
        {params.relate}/bin/Relate \
            --mode All \
            --haps {input.hapsout} \
            --sample {input.sampout} \
            -m {params.mu} \
            --coal {input.coa} \
            -o relate.mcmc.{params.mid}.temp \
            --map {input.genmap} \
            || true
        cp relate.mcmc.{params.mid}.temp.anc {output.anc}
        cp relate.mcmc.{params.mid}.temp.mut {output.mut}
        """

# resample the branch lengths
rule relate_branch:
    input:
        anc='{macro}/{micro}/{seed}/relate.mcmc.anc',
        mut='{macro}/{micro}/{seed}/relate.mcmc.mut',
        coa=macro+'/'+str(config['CHANGE']['RELATE']['COAL']),
    output:
        anc='{macro}/{micro}/{seed}/resample.mcmc.anc',
        mut='{macro}/{micro}/{seed}/resample.mcmc.mut',
        tmb='{macro}/{micro}/{seed}/resample.mcmc.timeb',
    params:
        relate=str(config['CHANGE']['FOLDERS']['SOFTWARE'])+'/'+str(config['CHANGE']['PROGRAMS']['RELATE']),
        mu=str(config['FIXED']['SIMULATE']['MU']),
        loc=str(config['FIXED']['SIMULATE']['LOC']),
        num=str(config['CHANGE']['RELATE']['NUMSAM']),
        mid='{macro}/{micro}/{seed}',
        mid2='{macro}.{micro}.{seed}',
    shell:
        """
        {params.relate}/scripts/SampleBranchLengths/SampleBranchLengths.sh \
            --input {params.mid}/relate.mcmc \
            --output {params.mid}/resample.mcmc.temp \
            -m {params.mu} \
            --coal {input.coa} \
            --format b \
            --num_samples {params.num} \
            --first_bp {params.loc} --last_bp {params.loc}
        cp {params.mid}/resample.mcmc.temp.anc {output.anc}
        cp {params.mid}/resample.mcmc.temp.mut {output.mut}
        cp {params.mid}/resample.mcmc.temp.timeb {output.tmb}
        rm -r relate.mcmc.{params.mid2}.temp/ || true
        rm relate.mcmc.{params.mid2}.temp.anc || true
        rm relate.mcmc.{params.mid2}.temp.mut || true
        rm -r {params.mid}/resample.mcmc.temp/ || true
        rm {params.mid}/resample.mcmc.temp.anc || true
        rm {params.mid}/resample.mcmc.temp.mut || true
        rm {params.mid}/resample.mcmc.temp.timeb || true
        rm {params.mid}/resample.mcmc.temp.dist || true
        """

# I use || to try to remove temporary files that could be large
# Relate doesn't play nice with working directory, hence the linux cp commands

# modify the inference file
rule clues:
    input:
        tmb='{macro}/{micro}/{seed}/resample.mcmc.timeb',
        frq='{macro}/{micro}/{seed}/slimulation.freq',
        coa=macro+'/'+str(config['CHANGE']['RELATE']['COAL']),
    output:
        anc='{macro}/{micro}/{seed}/clues.mcmc.epochs.npy',
        mut='{macro}/{micro}/{seed}/clues.mcmc.freqs.npy',
        pst='{macro}/{micro}/{seed}/clues.mcmc.post.npy',
        # modify clues inference.py to output this
        sel='{macro}/{micro}/{seed}/clues.mcmc.selcoef',
    params:
        freqscript=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        clues=str(config['CHANGE']['FOLDERS']['SOFTWARE'])+'/'+str(config['CHANGE']['PROGRAMS']['CLUES']),
        mu=str(config['FIXED']['SIMULATE']['MU']),
        loc=str(config['FIXED']['SIMULATE']['LOC']),
        thin=str(config['CHANGE']['CLUES']['THIN']),
        burn=str(config['CHANGE']['CLUES']['BURNIN']),
        mid='{macro}/{micro}/{seed}',
    shell:
        """
        daf=$(python {params.freqscript}/get-daf.py {input.frq})
        cd {params.clues}/
        python inference.py \
            --times {params.mid}/resample.mcmc \
            --coal {input.coa} \
            --out {params.mid}/clues.mcmc. \
            --popFreq $daf \
            --thin {params.thin} \
            --burnin {params.burn}
        rm {params.mid}/relate.mcmc.anc || true
        rm {params.mid}/relate.mcmc.mut || true
        rm {params.mid}/relate.mcmc.haps || true
        rm {params.mid}/relate.mcmc.sample || true
        rm {params.mid}/resample.mcmc.anc || true
        rm {params.mid}/resample.mcmc.mut || true
        """

# I have to put the clues/utils/ in the current directory

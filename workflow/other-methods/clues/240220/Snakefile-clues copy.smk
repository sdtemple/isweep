# isweep simulation studies - clues / relate
# Seth D. Temple, sdtemple@uw.edu
# May 23, 2023

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
sims['FILE'] = sims['FOLDER'] + '/isweep.inference.tsv'
sims['EXISTS'] = sims['FILE'].apply(os.path.isfile)
sims=sims[sims['EXISTS']]
sims = sims.set_index("SIMNAME", drop=False)

rule all:
	input:
		[f"{sim.FOLDER}/clues.selcoef" for sim in sims.itertuples()],
		[f"{sim.FOLDER}/relate.chr1.vcf.gz" for sim in sims.itertuples()],

##### relate ---> clues #####

# format coal file beforehand

# downsample, subset to region near causal mutation
rule relate_vcf:
    input:
        vcfin='{macro}/{micro}/{seed}/small.chr1.vcf.gz',
        subsample=macro+'/'+str(config['CHANGE']['RELATE']['SAMPLE']),
    output:
        vcfout='{macro}/{micro}/{seed}/relate.chr1.vcf.gz',
    params:
        loc=str(config['FIXED']['SIMULATE']['LOC']),
        pm=str(config['CHANGE']['RELATE']['BUFFER']),
    shell:
        """
        left=$(python -c "out = {params.loc} - {params.pm} ; print(out)")
        right=$(python -c "out = {params.loc} + {params.pm} ; print(out)")
        bcftools view {input.vcfin} \
            -r 1:${{left}}-${{right}} \
            -S {input.subsample}  \
            -Oz -o {output.vcfout}
        """

# set up input files from vcf
rule relate_input:
    input:
        vcfin='{macro}/{micro}/{seed}/relate.chr1.vcf.gz',
    output:
        hapsout='{macro}/{micro}/{seed}/relate.haps',
        sampout='{macro}/{micro}/{seed}/relate.sample',
    params:
        relate=str(config['CHANGE']['FOLDERS']['SOFTWARE'])+'/'+str(config['CHANGE']['PROGRAMS']['RELATE']),
    shell:
        """
        {params.relate}/bin/RelateFileFormats \
            --mode ConvertFromVcf \
            --haps {output.hapsout} \
            --sample {output.sampout} \
            --input {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/relate.chr1
        """

# run relate for first time
rule relate_run:
    input:
        hapsout='{macro}/{micro}/{seed}/relate.haps',
        sampout='{macro}/{micro}/{seed}/relate.sample',
        genmap='{macro}/uniform.map',
        coa=macro+'/'+str(config['CHANGE']['RELATE']['COAL']),
    output:
        anc='{macro}/{micro}/{seed}/relate.anc',
        mut='{macro}/{micro}/{seed}/relate.mut',
    params:
        relate=str(config['CHANGE']['FOLDERS']['SOFTWARE'])+'/'+str(config['CHANGE']['PROGRAMS']['RELATE']),
        mu=str(config['FIXED']['SIMULATE']['MU']),
    shell:
        """
        {params.relate}/bin/Relate \
            --mode All \
            --haps {input.hapsout} \
            --sample {input.sampout} \
            -m {params.mu} \
            --coal {input.coa} \
            -o relate.{wildcards.micro}.{wildcards.seed}.temp \
            --map {input.genmap} \
            || true
        cp relate.{wildcards.micro}.{wildcards.seed}.temp.anc {output.anc}
        cp relate.{wildcards.micro}.{wildcards.seed}.temp.mut {output.mut}
        """

# resample the branch lengths
rule relate_branch:
    input:
        anc='{macro}/{micro}/{seed}/relate.anc',
        mut='{macro}/{micro}/{seed}/relate.mut',
        coa=macro+'/'+str(config['CHANGE']['RELATE']['COAL']),
    output:
        anc='{macro}/{micro}/{seed}/resample.anc',
        mut='{macro}/{micro}/{seed}/resample.mut',
        tmb='{macro}/{micro}/{seed}/resample.timeb',
    params:
        relate=str(config['CHANGE']['FOLDERS']['SOFTWARE'])+'/'+str(config['CHANGE']['PROGRAMS']['RELATE']),
        mu=str(config['FIXED']['SIMULATE']['MU']),
        loc=str(config['FIXED']['SIMULATE']['LOC']),
        num=str(config['CHANGE']['RELATE']['NUMSAM'])
    shell:
        """
        {params.relate}/scripts/SampleBranchLengths/SampleBranchLengths.sh \
            --input {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/relate \
            --output {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/resample.temp \
            -m {params.mu} \
            --coal {input.coa} \
            --format b \
            --num_samples {params.num} \
            --first_bp {params.loc} --last_bp {params.loc}
        cp {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/resample.temp.anc {output.anc}
        cp {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/resample.temp.mut {output.mut}
        cp {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/resample.temp.timeb {output.tmb}
        rm -r relate.{wildcards.micro}.{wildcards.seed}.temp/ || true
        rm relate.{wildcards.micro}.{wildcards.seed}.temp.anc || true
        rm relate.{wildcards.micro}.{wildcards.seed}.temp.mut || true
        rm -r {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/resample.temp/ || true
        rm {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/resample.temp.anc || true
        rm {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/resample.temp.mut || true
        rm {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/resample.temp.timeb || true
        rm {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/resample.temp.dist || true
        """

# I use || to try to remove temporary files that could be large
# Relate doesn't play nice with working directory, hence the linux cp commands

# modify the inference file
rule clues:
    input:
        tmb='{macro}/{micro}/{seed}/resample.timeb',
        frq='{macro}/{micro}/{seed}/slimulation.freq',
        coa=macro+'/'+str(config['CHANGE']['RELATE']['COAL']),
    output:
        anc='{macro}/{micro}/{seed}/clues.epochs.npy',
        mut='{macro}/{micro}/{seed}/clues.freqs.npy',
        pst='{macro}/{micro}/{seed}/clues.post.npy',
        # modify clues inference.py to output this
        sel='{macro}/{micro}/{seed}/clues.selcoef',
    params:
        freqscript=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        clues=str(config['CHANGE']['FOLDERS']['SOFTWARE'])+'/'+str(config['CHANGE']['PROGRAMS']['CLUES']),
        mu=str(config['FIXED']['SIMULATE']['MU']),
        loc=str(config['FIXED']['SIMULATE']['LOC']),
        thin=str(config['CHANGE']['CLUES']['THIN']),
        burn=str(config['CHANGE']['CLUES']['BURNIN']),
    shell:
        """
        daf=$(python {params.freqscript}/get-daf.py {input.frq})
        python {params.clues}/inference.py \
            --times {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/resample \
            --coal {input.coa} \
            --out {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/clues \
            --popFreq $daf \
            --thin {params.thin} \
            --burnin {params.burn}
        rm {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/relate.anc || true
        rm {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/relate.mut || true
        rm {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/relate.haps || true
        rm {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/relate.sample || true
        rm {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/resample.anc || true
        rm {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/resample.mut || true
        """

# I have to put the clues/utils/ in the current directory

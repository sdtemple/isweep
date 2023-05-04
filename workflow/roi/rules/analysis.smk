# conduct analysis using isweep package
# seth temple, sdtemple@uw.edu
# may 3, 2023

# some input, string management, count sample size
subsamplefile=str(config['CHANGE']['ISWEEP']['SUBSAMPLE'])
macro=str(config['CHANGE']['FOLDERS']['STUDY'])
samplesize=0
with open(macro+'/'+subsamplefile,'r') as f:
    for line in f:
        samplesize+=1
samplesize=str(samplesize)

# rank polymorphisms in focus region
rule rank:
    input:
        short='{cohort}/{roi}/chr{chr}/center{center}/left{left}/right{right}/focus.ibd.gz',
        vcf='{cohort}/{roi}/chr{chr}/center{center}/left{left}/right{right}/chr.focused.vcf.gz',
    output:
        fileout='{cohort}/{roi}/chr{chr}/center{center}/left{left}/right{right}/isweep.ranks.tsv.gz',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        diameter=str(config['FIXED']['ISWEEP']['DIAMETER']),
        atleast=str(config['FIXED']['ISWEEP']['ATLEAST']),
        q1=str(config['FIXED']['ISWEEP']['MINAAF']),
        q2=str(config['FIXED']['ISWEEP']['MAXAAF']),
        rulesigma=str(config['FIXED']['ISWEEP']['RULESIGMA']),
    shell:
        """
        python {params.scripts}/rank-isweep.py \
            {input.short} \
            {input.vcf} \
            {output.fileout} \
            {params.diameter} \
            {params.atleast} \
            {params.q1} \
            {params.q2} \
            {params.rulesigma}
        """
        # 'python {config[CHANGE][FOLDERS][TERMINALSCRIPTS]}/rank-isweep.py {input.short} {input.vcf} {output.fileout} {params.diameter} {params.atleast} {params.q1} {params.q2} {params.rulesigma}'

# some input, string management
macro=str(config['CHANGE']['FOLDERS']['STUDY'])
micro=str(config['CHANGE']["ISWEEP"]["ROI"])
sims = pd.read_csv(macro+'/'+micro, sep='\t', header=0)
sims['FOLDER'] = [(macro +'/'+str(sims.loc[j].NAME)+'/chr'+str(sims.iloc[j].CHROM)+'/center'+str(sims.loc[j].BPCENTER)+'/left'+str(sims.loc[j].BPLEFT)+'/right'+str(sims.loc[j].BPRIGHT)).strip() for j in range(J)]

# extend Ne(t)
rule extendNe:
    input:
        shortNe=macro+'/'+str(config['CHANGE']['ISWEEP']['NE']),
    output:
        longNe=macro+'/extended.ne',
    shell:
        """
        python -c 'from isweep import extend_Ne; extend_Ne("{input.shortNe}",1000,"{output.longNe}")'
        """

# make isweep inferences
rule infer:
    input:
        long='{cohort}/{roi}/chr{chr}/center{center}/left{left}/right{right}/mom.ibd.gz',
        ranks='{cohort}/{roi}/chr{chr}/center{center}/left{left}/right{right}/isweep.ranks.tsv.gz',
        effdemo=macro+'/extended.ne',
    output:
        fileout='{cohort}/{roi}/chr{chr}/center{center}/left{left}/right{right}/isweep.inference.tsv',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        nboot=str(config['FIXED']['ISWEEP']['NBOOT']),
        cm=str(config['FIXED']['ISWEEP']['MOMCUTOFF']),
        n=samplesize,
        ploidy=str(config['FIXED']['HAPIBD']['PLOIDY']),
    shell:
        """
        python {params.scripts}/infer-isweep.py \
            {input.long} \
            {input.ranks} \
            {output.fileout} \
            {params.nboot} \
            {params.cm} \
            {params.n} \
            {input.effdemo} \
            {params.ploidy}
        """
        # 'python {config[CHANGE][FOLDERS][TERMINALSCRIPTS]}/infer-isweep.py {input.long} {input.ranks} {output.fileout} {params.nboot} {params.cm} {params.n} {params.effdemo} {params.ploidy}'

rule combine:
    input:
        sweeps=[f"{sim.FOLDER}/isweep.inference.tsv" for sim in sims.itertuples()],
        roitsv=macro+'/'+micro,
    output:
        doneso=macro+'/doneso.txt',
    shell:
        'touch {output.doneso}'

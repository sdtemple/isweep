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
        q1=str(config['FIXED']['ISWEEP']['MINMAXAAF']),
        rulesigma=str(config['FIXED']['ISWEEP']['RULESIGMA']),
    shell:
        """
        python {params.scripts}/rank-isweep.py \
            {input.short} \
            {input.vcf} \
            {output.fileout} \
            {params.diameter} \
            {params.q1} \
            {params.rulesigma}
        """

# some input, string management
macro=str(config['CHANGE']['FOLDERS']['STUDY'])
micro=str(config['CHANGE']["ISWEEP"]["ROI"])
sims = pd.read_csv(macro+'/'+micro, sep='\t', header=0)
sims['FOLDER'] = [(macro +'/'+str(sims.loc[j].NAME)+'/chr'+str(sims.iloc[j].CHROM)+'/center'+str(sims.loc[j].BPCENTER)+'/left'+str(sims.loc[j].BPLEFTCENTER)+'/right'+str(sims.loc[j].BPRIGHTCENTER)).strip() for j in range(J)]

# # extend Ne(t)
# rule extendNe:
#     input:
#         shortNe=macro+'/'+str(config['CHANGE']['ISWEEP']['NE']),
#     output:
#         longNe=macro+'/extended.ne',
#     shell:
#         """
#         python -c 'from isweep import extend_Ne; extend_Ne("{input.shortNe}",1000,"{output.longNe}")'
#         """

# # make isweep inferences
# rule infer:
#     input:
#         long='{cohort}/{roi}/chr{chr}/center{center}/left{left}/right{right}/mom.ibd.gz',
#         ranks='{cohort}/{roi}/chr{chr}/center{center}/left{left}/right{right}/isweep.ranks.tsv.gz',
#         effdemo=macro+'/extended.ne',
#     output:
#         fileout='{cohort}/{roi}/chr{chr}/center{center}/left{left}/right{right}/isweep.inference.tsv',
#     params:
#         scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
#         nboot=str(config['FIXED']['ISWEEP']['NBOOT']),
#         cm=str(config['FIXED']['ISWEEP']['MOMCUTOFF']),
#         n=str(samplesize),
#         ploidy=str(config['FIXED']['HAPIBD']['PLOIDY']),
#         quant=str(config['FIXED']['ISWEEP']['QUANT']),
#     shell:
#         """
#         python {params.scripts}/infer-isweep.py \
#             {input.long} \
#             {input.ranks} \
#             {output.fileout} \
#             {params.nboot} \
#             {params.cm} \
#             {params.n} \
#             {input.effdemo} \
#             {params.ploidy} \
#             {params.quant}
#         """

rule combine:
    input:
        # sweeps=[f"{sim.FOLDER}/isweep.inference.tsv" for sim in sims.itertuples()],
        ranks=[f"{sim.FOLDER}/isweep.ranks.tsv.gz" for sim in sims.itertuples()],
	output:
		yamlout=macro+'/arguments.roi.yaml',
	params:
		yamlpar=macro+'/'+yaml,
	shell:
		'cp {params.yamlpar} {output.yamlout}'

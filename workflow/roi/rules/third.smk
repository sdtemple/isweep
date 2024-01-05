### estimating selection coefficients
### estimating allele frequencies

# some inputs, string managements, count sample size
subsamplefile=str(config['CHANGE']['ISWEEP']['SUBSAMPLE'])
cohort=str(config['CHANGE']['FOLDERS']['STUDY'])
samplesize=0
with open(cohort+'/'+subsamplefile,'r') as f:
    for line in f:
        samplesize+=1
ploidy=2
# ploidy=int(float(str(config['FIXED']['HAPIBD']['PLOIDY'])))
samplesize_ploidy=samplesize*ploidy

# extend Ne(t)
rule third_Ne:
    input:
        shortNe='{cohort}/'+str(config['CHANGE']['ISWEEP']['NE']),
    output:
        longNe='{cohort}/extended.ne',
    params:
        lastne=str(config['FIXED']['ISWEEP']['LASTNE']),
    shell:
        """
        python -c "from isweep import extend_Ne; extend_Ne('{input.shortNe}',float('{params.lastne}'),'{output.longNe}')"
        """

##### haplotype analysis #####

rule third_hap:
    input:
        rankin='{cohort}/{hit}/second.ranks.tsv.gz',
        outlied='{cohort}/{hit}/second.outliers.txt',
    output:
        lociout='{cohort}/{hit}/third.best.hap.txt',
        happng='{cohort}/{hit}/third.hap.png',
        snppng='{cohort}/{hit}/third.snp.png',
    params:
        windowsize=str(config['FIXED']['ISWEEP']['HAPSIZE']),
        windowstep=str(config['FIXED']['ISWEEP']['HAPSTEP']),
        freqsize=str(config['FIXED']['ISWEEP']['FREQSIZE']),
        freqstep=str(config['FIXED']['ISWEEP']['FREQSTEP']),
        numsnp=str(config['FIXED']['ISWEEP']['NUMSNP']),
        lowbnd=str(config['FIXED']['ISWEEP']['MINAAF']),
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
    shell:
        """
        python {params.scripts}/haplotypes.py \
            {input.rankin} \
            {wildcards.cohort}/{wildcards.hit} \
            0 \
            1 \
            -1 \
            {params.freqsize} \
            {params.freqstep} \
            {params.windowsize} \
            {params.windowstep} \
            {params.numsnp} \
            {params.lowbnd}
        """


rule third_hap_ibd:
    input:
        best='{cohort}/{hit}/third.best.hap.txt',
        locus='{cohort}/{hit}/locus.txt',
    output:
        ibd='{cohort}/{hit}/third.hap.ibd.gz',
    params:
        ibdfolder='{cohort}/ibdsegs/ibdends/mle',
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['FILTER']),
        script=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS'])+'/lines.py',
    shell:
        """
        chr=$(python {params.script} {input.locus} 2 2)
        thecenter=$(python {params.script} {input.best} 1 2)
        ibd={params.ibdfolder}/chr${{chr}}.ibd.gz
        zcat $ibd | \
            java -Xmx{config[CHANGE][ISWEEP][XMXMEM]}g -jar {params.soft}/{params.prog} \
            "I" 6 0.00 ${{thecenter}} | \
            java -Xmx{config[CHANGE][ISWEEP][XMXMEM]}g -jar {params.soft}/{params.prog} \
            "I" 7 ${{thecenter}} 10000000000 | \
            gzip > {output.ibd}
        """


rule third_hap_infer:
    input:
        long='{cohort}/{hit}/third.hap.ibd.gz',
        freq='{cohort}/{hit}/third.best.hap.txt',
        longNe='{cohort}/extended.ne',
    output:
        fileout='{cohort}/{hit}/results.hap.tsv',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        nboot=str(config['FIXED']['ISWEEP']['NBOOT']),
        cm=str(config['FIXED']['ISWEEP']['MLECUTOFF']),
        n=str(samplesize),
        ploidy=str(2),
        # ploidy=str(config['FIXED']['HAPIBD']['PLOIDY'])
    shell:
        """
        ibdest=$(zcat {input.long} | wc -l)
        freqest=$(python {params.scripts}/lines.py {input.freq} 2 2)
        python {params.scripts}/estimate.py \
            {output.fileout} \
            ${{ibdest}} \
            ${{freqest}} \
            {params.nboot} \
            {params.cm} \
            {params.n} \
            {input.longNe} \
            {params.ploidy}
        """

##### snp analysis #####

rule third_snp:
    input:
        rankin='{cohort}/{hit}/second.ranks.tsv.gz',
        outlied='{cohort}/{hit}/second.outliers.txt',
    output:
        lociout='{cohort}/{hit}/third.best.snp.txt',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        lowbnd=str(config['FIXED']['ISWEEP']['MINAAF']),
    shell:
        """
        python {params.scripts}/snp.py \
            {input.rankin} \
            {output.lociout} \
            {params.lowbnd}
        """

rule third_snp_ibd:
    input:
        best='{cohort}/{hit}/third.best.snp.txt',
        locus='{cohort}/{hit}/locus.txt',
    output:
        ibd='{cohort}/{hit}/third.snp.ibd.gz',
    params:
        ibdfolder='{cohort}/ibdsegs/ibdends/mle',
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['FILTER']),
        script=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS'])+'/lines.py',
    shell:
        """
        chr=$(python {params.script} {input.locus} 2 2)
        thecenter=$(python {params.script} {input.best} 1 2)
        ibd={params.ibdfolder}/chr${{chr}}.ibd.gz
        zcat $ibd | \
            java -Xmx{config[CHANGE][ISWEEP][XMXMEM]}g -jar {params.soft}/{params.prog} \
            "I" 6 0.00 ${{thecenter}} | \
            java -Xmx{config[CHANGE][ISWEEP][XMXMEM]}g -jar {params.soft}/{params.prog} \
            "I" 7 ${{thecenter}} 10000000000 | \
            gzip > {output.ibd}
        """


rule third_snp_infer:
    input:
        long='{cohort}/{hit}/third.snp.ibd.gz',
        freq='{cohort}/{hit}/third.best.snp.txt',
        longNe='{cohort}/extended.ne',
    output:
        fileout='{cohort}/{hit}/results.snp.tsv',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        nboot=str(config['FIXED']['ISWEEP']['NBOOT']),
        cm=str(config['FIXED']['ISWEEP']['MLECUTOFF']),
        n=str(samplesize),
        ploidy=str(2),
        # ploidy=str(config['FIXED']['HAPIBD']['PLOIDY'])
    shell:
        """
        ibdest=$(zcat {input.long} | wc -l)
        freqest=$(python {params.scripts}/lines.py {input.freq} 2 2)
        python {params.scripts}/estimate.py \
            {output.fileout} \
            ${{ibdest}} \
            ${{freqest}} \
            {params.nboot} \
            {params.cm} \
            {params.n} \
            {input.longNe} \
            {params.ploidy}
        """

### write entropy ###

rule entropy:
	input:
		filein='{cohort}/{hit}/outlier1.txt',
	output:
		fileout='{cohort}/{hit}/ibd.entropy.tsv',
	params:
		scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
		samplesizep=str(samplesize_ploidy),
	shell:
		"""
		python {params.scripts}/ibd-entropy.py \
			{wildcards.cohort}/{wildcards.hit} \
			{output.fileout} \
			{params.samplesizep}
		"""

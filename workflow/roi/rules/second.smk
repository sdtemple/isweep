wildcard_constraints:
	NAME = '\w+',

# some inputs, string managements, count sample size
cohort=str(config['CHANGE']['FOLDERS']['STUDY'])
ploidy=int(float(config['FIXED']['HAPIBD']['PLOIDY']))
maf=float(config['FIXED']['ISWEEP']['MINAAF2'])

# subset vcf to region of interest
rule second_region: # focus vcf on region of interest
    input:
        locus='{cohort}/{hit}/locus.txt',
        focus='{cohort}/{hit}/first.pos.txt',
        subsample=cohort+"/subsample.txt",
        excludesamples=cohort+"/excludesamples.txt",
    output:
        subvcf='{cohort}/{hit}/second.focused.vcf.gz',
    params:
        qmaf=maf,
        chrpre=str(config['CHANGE']['ISWEEP']['CHRPRE']),
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        vcfs=str(config['CHANGE']['EXISTING']['VCFS']),
        pm=str(config['FIXED']['ISWEEP']['PM']),
    resources:
        mem_gb=10
    shell: # if chromosome is huge (greater than 10000 Mb), may need to modify the third pipe
        """
        chr=$(python {params.scripts}/lines.py {input.locus} 2 2)
        center=$(python {params.scripts}/lines.py {input.focus} 1 2)
        left=$(python -c "out = $center - {params.pm} ; print(out)")
        right=$(python -c "out = $center + {params.pm} ; print(out)")
        vcf={params.vcfs}/chr${{chr}}.vcf.gz
        tabix -fp vcf ${{vcf}}
        bcftools view ${{vcf}} -r {params.chrpre}${{chr}}:${{left}}-${{right}} -Ob | \
            bcftools view -S {input.subsample} -Ob | \
            bcftools view -q {params.qmaf}:nonmajor -Oz -o {output.subvcf}
        """

### filter ibd file ###

rule second_filt:
    input:
        ibd='{cohort}/{hit}/first.ibd.gz',
        locus='{cohort}/{hit}/first.pos.txt',
    output:
        ibd='{cohort}/{hit}/second.filt.ibd.gz',
    params:
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['FILTER']),
        script=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS'])+'/lines.py',
    resources:
        mem_gb='{config[CHANGE][ISWEEP][XMXMEM]}'
    shell:
        """
        thecenter=$(python {params.script} {input.locus} 1 2)
        zcat {input.ibd} | \
            java -Xmx{config[CHANGE][ISWEEP][XMXMEM]}g -jar {params.soft}/{params.prog} \
            "I" 6 0.00 $thecenter | \
            java -Xmx{config[CHANGE][ISWEEP][XMXMEM]}g -jar {params.soft}/{params.prog} \
            "I" 7 $thecenter 10000000000 | \
            gzip > {output.ibd}
        """

### rank snps ###

rule second_rank:
    input:
        short='{cohort}/{hit}/second.filt.ibd.gz',
        vcf='{cohort}/{hit}/second.focused.vcf.gz',
    output:
        fileout='{cohort}/{hit}/second.ranks.tsv.gz',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        diameter=str(config['FIXED']['ISWEEP']['DIAMETER']),
        q1=str(config['FIXED']['ISWEEP']['MINAAF2']),
        rulesigma=str(config['FIXED']['ISWEEP']['RULESIGMA']),
    shell:
        """
        python {params.scripts}/rank.py \
            {input.short} \
            {input.vcf} \
            {output.fileout} \
            {params.diameter} \
            {params.q1} \
            {params.rulesigma}
        """

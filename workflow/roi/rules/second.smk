### finding ibd groups

# some inputs, string managements, count sample size
cohort=str(config['CHANGE']['FOLDERS']['STUDY'])
ploidy=2
# ploidy=int(float(config['FIXED']['HAPIBD']['PLOIDY']))
maf=float(config['FIXED']['ISWEEP']['MINAAF'])

# subset vcf to region of interest
rule second_region: # focus vcf on region of interest
    input:
        locus='{cohort}/{hit}/locus.txt',
        focus='{cohort}/{hit}/first.pos.txt',
        vcf='{cohort}/{hit}/first.focused.vcf.gz',
        subsample=cohort+"/subsample.txt",
    output:
        subvcf='{cohort}/{hit}/second.focused.vcf.gz',
    params:
        qmaf=maf,
        chrpre=str(config['CHANGE']['ISWEEP']['CHRPRE']),
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        pm=str(config['FIXED']['ISWEEP']['PM']),
    shell: # if chromosome is huge (greater than 10000 Mb), may need to modify the third pipe
        """
        chr=$(python {params.scripts}/lines.py {input.locus} 2 2)
        center=$(python {params.scripts}/lines.py {input.focus} 1 2)
        left=$(python -c "out = $center - {params.pm} ; print(max(out,2))")
        right=$(python -c "out = $center + {params.pm} ; print(out)")
        tabix -fp vcf {input.vcf}
        bcftools view {input.vcf} -r {params.chrpre}${{chr}}:${{left}}-${{right}} -Ob | \
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
        q1=str(config['FIXED']['ISWEEP']['MINAAF']),
        rulesigma=str(config['FIXED']['ISWEEP']['GROUPCUTOFF']),
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

### write outliers ###

rule second_outlier:
    input:
        short='{cohort}/{hit}/second.filt.ibd.gz',
    output:
        fileout='{cohort}/{hit}/second.outliers.txt',
        out1='{cohort}/{hit}/outlier1.txt',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        diameter=str(config['FIXED']['ISWEEP']['DIAMETER']),
        rulesigma=str(config['FIXED']['ISWEEP']['GROUPCUTOFF']),
    shell:
        """
        python {params.scripts}/outliers.py \
            {input.short} \
            {wildcards.cohort}/{wildcards.hit} \
            {params.diameter} \
            {params.rulesigma}
        touch {output.fileout}
        touch {output.out1}
        """

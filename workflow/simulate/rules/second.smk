### refining location of analysis

# inputs
n=int(float(config['CHANGE']['SIMULATE']['SAMPSIZE']))
ploidy=int(float(config['FIXED']['SIMULATE']['PLOIDY']))
maf3=float(config['FIXED']['HAPIBD']['MINMAF'])
mac3=int(ploidy*n*maf3)

### filter vcf ###

rule second_region:
    input:
        vcfin='{macro}/{micro}/{seed}/large.chr1.vcf.gz',
        locus='{macro}/{micro}/{seed}/first.pos.txt',
    output:
        vcfout='{macro}/{micro}/{seed}/second.chr1.vcf.gz',
    params:
        folder='{macro}/{micro}/{seed}',
        pm=str(config['FIXED']['SIMULATE']['BUFFER']),
        script=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS'])+'/lines.py',
    shell: # to bgz and back is being consertative
        """
        thecenter=$(python ../../scripts/lines.py {input.locus} 1 2)
        gunzip -c {input.vcfin} | bgzip  > {params.folder}/chrtemp.vcf.bgz
        tabix -fp vcf {params.folder}/chrtemp.vcf.bgz
        left=$(python -c "out = $thecenter - {params.pm} ; print(out)")
        right=$(python -c "out = $thecenter + {params.pm} ; print(out)")
        bcftools view {params.folder}/chrtemp.vcf.bgz \
            -r 1:${{left}}-${{right}} \
            -Oz -o {output.vcfout}
        rm {params.folder}/chrtemp.vcf.bgz
        """

### call hap-ibd ###

rule second_hapibd:
    input:
        vcf='{macro}/{micro}/{seed}/second.chr1.vcf.gz',
        map=str(config["CHANGE"]["FOLDERS"]["MACRO"]) + '/uniform.map',
    params:
        minmac=str(mac3),
        out='{macro}/{micro}/{seed}/second.chr1',
        minsee=str(config['FIXED']['HAPIBD']['MINSEED']),
        minext=str(config['FIXED']['HAPIBD']['MINEXT']),
        minout=str(config['FIXED']['HAPIBD']['MINOUT']),
    output:
        ibd='{macro}/{micro}/{seed}/second.chr1.ibd.gz',
        hbd='{macro}/{micro}/{seed}/second.chr1.hbd.gz',
        log='{macro}/{micro}/{seed}/second.chr1.log',
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}',
    shell:
        """
        java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar ../../software/hap-ibd.jar \
            gt={input.vcf} \
            map={input.map} \
            out={params.out} \
            min-seed={params.minsee} \
            min-extend={params.minext} \
            min-output={params.minout} \
            min-mac={params.minmac}
        """

### filter ibd file ###

rule second_filt:
    input:
        ibd='{macro}/{micro}/{seed}/second.chr1.ibd.gz',
        locus='{macro}/{micro}/{seed}/first.pos.txt',
    output:
        ibd='{macro}/{micro}/{seed}/second.filt.chr1.ibd.gz',
    shell:
        """
        thecenter=$(python ../../scripts/lines.py {input.locus} 1 2)
        python ../../scripts/filter-lines.py \
            --file_input {input.ibd} \
            --file_output {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/intermediate.ibd.gz \
            --column_index 6 \
            --upper_bound $thecenter \
            --complement 0
        python ../../scripts/filter-lines.py \
            --file_input {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/intermediate.ibd.gz \
            --file_output {output.ibd} \
            --column_index 7 \
            --lower_bound $thecenter \
            --upper_bound 10000000000 \
            --complement 0
        rm {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/intermediate.ibd.gz
        """

### rank snps ###

rule second_rank:
    input:
        short='{macro}/{micro}/{seed}/second.filt.chr1.ibd.gz',
        vcf='{macro}/{micro}/{seed}/second.chr1.vcf.gz',
    output:
        fileout='{macro}/{micro}/{seed}/second.ranks.tsv.gz',
    params:
        diameter=str(config['FIXED']['ISWEEP']['DIAMETER']),
        q1=str(config['FIXED']['ISWEEP']['MINAAF']),
        rulesigma=str(config['FIXED']['ISWEEP']['GROUPCUTOFF']),
    shell:
        """
        python ../../scripts/rank.py \
            --ibd_file {input.short} \
            --vcf {input.vcf} \
            --file_out {output.fileout} \
            --graph_diameter {params.diameter} \
            --low_freq {params.q1} \
            --group_cutoff {params.rulesigma} \
        """
        
### write outliers ###

rule second_outlier:
    input:
        short='{macro}/{micro}/{seed}/second.filt.chr1.ibd.gz',
    output:
        fileout='{macro}/{micro}/{seed}/second.outliers.txt',
    params:
        diameter=str(config['FIXED']['ISWEEP']['DIAMETER']),
        rulesigma=str(config['FIXED']['ISWEEP']['GROUPCUTOFF']),
    shell:
        """
        python ../../scripts/outliers.py \
            --file_input {input.short} \
            --folder_output {wildcards.macro}/{wildcards.micro}/{wildcards.seed} \
            --graph_diameter {params.diameter} \
            --group_cutoff {params.rulesigma}
        touch {output.fileout}
        """

rule gini_impurity:
	input:
		filein='{macro}/{micro}/{seed}/outlier1.txt',
	output:
		fileout='{macro}/{micro}/{seed}/ibd.gini.tsv',
	params:
		samplesizep=str(samplesize_ploidy),
	shell:
		"""
		python ../../scripts/ibd-gini-entropy.py \
			--folder {wildcards.cohort}/{wildcards.hit} \
			--file_out {output.fileout} \
			--sample_size {params.samplesizep}
		"""
### refining location of analysis

# inputs
n=int(float(config['CHANGE']['SIMULATE']['SAMPSIZE']))
ploidy=int(float(config['FIXED']['SIMULATE']['PLOIDY']))
maf3=float(config['FIXED']['HAPIBD']['MINMAF'])
mac3=int(ploidy*n*maf3)
samplesize_ploidy=n*ploidy

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
    shell: # to bgz and back is being consertative
        """
        thecenter=$(python ../../scripts/utilities/lines.py {input.locus} 1 2)
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
        thecenter=$(python ../../scripts/utilities/lines.py {input.locus} 1 2)
        python ../../scripts/utilities/filter-lines.py \
            --input_file {input.ibd} \
            --output_file {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/intermediate.ibd.gz \
            --column_index 6 \
            --upper_bound $thecenter \
            --complement 0
        python ../../scripts/utilities/filter-lines.py \
            --input_file {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/intermediate.ibd.gz \
            --output_file {output.ibd} \
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
        python ../../scripts/model/rank.py \
            --input_ibd_file {input.short} \
            --input_vcf_file {input.vcf} \
            --output_file {output.fileout} \
            --graph_diameter {params.diameter} \
            --lowest_freq {params.q1} \
            --group_cutoff {params.rulesigma} \
        """
        
### write outliers ###

rule second_outlier:
    input:
        short='{macro}/{micro}/{seed}/second.filt.chr1.ibd.gz',
    output:
        fileout='{macro}/{micro}/{seed}/second.outliers.txt',
        out1='{macro}/{micro}/{seed}/outlier1.txt',
    params:
        diameter=str(config['FIXED']['ISWEEP']['DIAMETER']),
        rulesigma=str(config['FIXED']['ISWEEP']['GROUPCUTOFF']),
    shell:
        """
        python ../../scripts/model/outliers.py \
            --input_ibd_file {input.short} \
            --output_folder {wildcards.macro}/{wildcards.micro}/{wildcards.seed} \
            --graph_diameter {params.diameter} \
            --group_cutoff {params.rulesigma}
        touch {output.fileout}
        touch {output.out1}
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
		python ../../scripts/model/ibd-gini-entropy.py \
			--input_folder {wildcards.macro}/{wildcards.micro}/{wildcards.seed} \
			--output_file {output.fileout} \
			--sample_size {params.samplesizep}
		"""
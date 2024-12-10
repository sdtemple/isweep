### initializing regions of interest
### some localizing

# inputs
n=int(float(config['CHANGE']['SIMULATE']['SAMPSIZE']))
ploidy=int(float(config['FIXED']['SIMULATE']['PLOIDY']))
maf3=float(config['FIXED']['HAPIBD']['MINMAF'])
mac3=int(ploidy*n*maf3)

### modal location

rule first_mode:
    input:
        ibd='{macro}/{micro}/{seed}/scan.chr1.windowed.tsv.gz',
    output:
        ibd='{macro}/{micro}/{seed}/first.mode.txt',
    shell:
        """
        python scripts/where-mode-ibd.py \
            {input.ibd} \
            {output.ibd} \
        """

### filter vcf ###

rule first_region:
    input:
        vcfin='{macro}/{micro}/{seed}/large.chr1.vcf.gz',
        ibdwin='{macro}/{micro}/{seed}/first.mode.txt',
    output:
        vcfout='{macro}/{micro}/{seed}/first.chr1.vcf.gz',
    params:
        folder='{macro}/{micro}/{seed}',
        pm=str(config['FIXED']['SIMULATE']['BUFFER']),
    shell: # to bgz and back is being consertative
        """
        thecenter=$(python ../../scripts/utilities/lines.py {input.ibdwin} 1 2)
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

rule first_hapibd:
    input:
        vcf='{macro}/{micro}/{seed}/first.chr1.vcf.gz',
        map=str(config["CHANGE"]["FOLDERS"]["MACRO"]) + '/uniform.map',
    params:
        minmac=str(mac3),
        out='{macro}/{micro}/{seed}/first.chr1',
        minsee=str(config['FIXED']['HAPIBD']['MINSEED']),
        minext=str(config['FIXED']['HAPIBD']['MINEXT']),
        minout=str(config['FIXED']['HAPIBD']['MINOUT']),
    output:
        ibd='{macro}/{micro}/{seed}/first.chr1.ibd.gz',
        hbd='{macro}/{micro}/{seed}/first.chr1.hbd.gz',
        log='{macro}/{micro}/{seed}/first.chr1.log',
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

rule first_filt:
    input:
        ibd='{macro}/{micro}/{seed}/first.chr1.ibd.gz',
        locus='{macro}/{micro}/{seed}/first.mode.txt',
    output:
        ibd='{macro}/{micro}/{seed}/first.filt.chr1.ibd.gz',
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

rule first_rank:
    input:
        short='{macro}/{micro}/{seed}/first.filt.chr1.ibd.gz',
        vcf='{macro}/{micro}/{seed}/first.chr1.vcf.gz',
    output:
        fileout='{macro}/{micro}/{seed}/first.ranks.tsv.gz',
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

rule first_score:
    input:
        snps='{macro}/{micro}/{seed}/first.ranks.tsv.gz',
    output:
        loci='{macro}/{micro}/{seed}/first.pos.txt',
    params:
        windowsize=str(config['FIXED']['ISWEEP']['WINSIZE']),
        windowstep=str(config['FIXED']['ISWEEP']['WINSTEP']),
        qrng=str(config['FIXED']['ISWEEP']['QRANGE']),
        maxspace=str(config['FIXED']['ISWEEP']['MAXSPACING']),
        folderout='{macro}/{micro}/{seed}',
    shell:
        """
        python ../../scripts/model/site.py \
            --input_snp_file {input.snps} \
            --output_folder {params.folderout} \
            --window_index 0 \
            --freq_index 1 \
            --window_size {params.windowsize} \
            --window_step {params.windowstep} \
            --quantile_sum {params.qrng} \
            --max_spacing {params.maxspace} \
        """

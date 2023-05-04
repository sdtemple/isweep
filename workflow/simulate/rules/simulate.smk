wildcard_constraints:
    SIMNAME = '\w+',

rule uniform_map:
    input:
    output:
        mapout=str(config["CHANGE"]["FOLDERS"]["MACRO"]) + '/uniform.map',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['SNAKESCRIPTS']),
    script:
        '{params.scripts}/uniformMap.py'

rule forward_Ne:
    input:
        str(config["CHANGE"]["FOLDERS"]["MACRO"]) + '/uniform.map',
    output:
        str(config["CHANGE"]["FOLDERS"]["MACRO"]) + '/forward.ne',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['SNAKESCRIPTS']),
    script:
        '{params.scripts}/forwardNe.py'

rule slim_script:
    input:
        trueNe = str(config["CHANGE"]["FOLDERS"]["MACRO"]) + '/forward.ne',
        map = str(config["CHANGE"]["FOLDERS"]["MACRO"]) + '/uniform.map',
    output:
        [f"{sim.FOLDER}/slimulation.slim" for sim in sims.itertuples()],
    params:
        scripts=str(config['CHANGE']['FOLDERS']['SNAKESCRIPTS']),
    script:
        '{params.scripts}/writeSlimDemography.py'

### generate tree, vcf data ###

# slim forward
rule slim:
    input:
        "{macro}/{micro}/{seed}/slimulation.slim",
    output:
        trees="{macro}/{micro}/{seed}/slimulation.trees",
        freq="{macro}/{micro}/{seed}/slimulation.freq",
    params:
        scripts=str(config['CHANGE']['FOLDERS']['SNAKESCRIPTS']),
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['SLiM']),
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}',
    shell:
        '{params.soft}/{params.prog} {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/slimulation.slim'

# msprime backward
rule msprime:
    input:
        trees = "{macro}/{micro}/{seed}/slimulation.trees",
    output:
        bcf="{macro}/{micro}/{seed}/slimulation.bcf.gz",
    params:
        scripts=str(config['CHANGE']['FOLDERS']['SNAKESCRIPTS']),
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}',
    script:
        '{params.scripts}/treeVCF.py'

# bcf, vcf magic
rule vcf:
# rule tabix_bcf:
    input:
        bcf='{macro}/{micro}/{seed}/slimulation.bcf.gz',
    output:
        vcf='{macro}/{micro}/{seed}/slimulation.vcf.gz',
        csi='{macro}/{micro}/{seed}/slimulation.bcf.gz.csi',
    shell:
        """
        tabix {input.bcf}
        bcftools view {input.bcf} -Oz -m2 -M2 -o {output.vcf}
        """

rule genotyping_error:
    input:
        vcf='{macro}/{micro}/{seed}/slimulation.vcf.gz',
    output:
        out='{macro}/{micro}/{seed}/large.chr1.vcf.gz',
    params:
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['GTERR']),
        gter=str(config['FIXED']['SIMULATE']['GTERR']),
    shell:
        """
        zcat {input.vcf} | \
            java -jar {params.soft}/{params.prog} {params.gter} | \
            gzip > {output.out}
        """
		# 'zcat {input.vcf} | java -jar {params.soft}/{params.prog} {params.gter} | gzip > {output.out}'

### true labels ###

# work on location
rule causal_vcf:
    input:
        bcf='{macro}/{micro}/{seed}/slimulation.bcf.gz',
        csi='{macro}/{micro}/{seed}/slimulation.bcf.gz.csi',
    output:
        vcf='{macro}/{micro}/{seed}/causal.chr1.vcf.gz'
    params:
        loc=str(config['FIXED']['SIMULATE']['LOC'])
    shell:
        'bcftools view {input.bcf} -Oz -m2 -M2 -t 1:{params.loc} -o {output.vcf}'

rule true_labels:
    input:
        vcf='{macro}/{micro}/{seed}/causal.chr1.vcf.gz'
    output:
        labels1='{macro}/{micro}/{seed}/slimulation.true1.labels.gz',
        labels0='{macro}/{micro}/{seed}/slimulation.true0.labels.gz',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['SNAKESCRIPTS']),
    script:
        '{params.scripts}/true-labels.py'

### zip tree sequence ###

rule tszip:
    input:
        trees='{macro}/{micro}/{seed}/slimulation.trees',
        yaml='{macro}/{micro}/{seed}/slimulation.bcf.gz',
    output:
        tsz='{macro}/{micro}/{seed}/slimulation.trees.tsz'
    shell:
        'tszip {input.trees}'

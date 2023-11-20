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
        gter=str(config['CHANGE']['SIMULATE']['GTERR']),
    shell:
        """
        zcat {input.vcf} | \
            java -jar {params.soft}/{params.prog} {params.gter} | \
            gzip > {output.out}
        """

### zip tree sequence ###

rule tszip:
    input:
        trees='{macro}/{micro}/{seed}/slimulation.trees',
        yaml='{macro}/{micro}/{seed}/slimulation.bcf.gz',
    output:
        tsz='{macro}/{micro}/{seed}/slimulation.trees.tsz'
    shell:
        'tszip {input.trees}'

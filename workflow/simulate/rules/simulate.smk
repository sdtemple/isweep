### simulate data

# inputs
ploidy=int(float(config['FIXED']['SIMULATE']['PLOIDY']))

rule uniform_map:
    input:
    output:
        mapout=str(config["CHANGE"]["FOLDERS"]["MACRO"]) + '/uniform.map',
    params:
        rho=str(config['CHANGE']['SIMULATE']['RHO']),
        L=str(config['CHANGE']['SIMULATE']['CMLEN']),
    shell:
        """
        python scripts/uniformMap.py \
            {output.mapout} \
            {params.rho} \
            {params.L}
        """

rule forward_Ne:
    input:
        nefile=str(config['CHANGE']['SIMULATE']['tNe']),
    output:
        fileout=str(config["CHANGE"]["FOLDERS"]["MACRO"]) + '/forward.ne',
    shell:
        """
        python scripts/forwardNe.py \
            {input.nefile} \
            {output.fileout}
        """

rule slim_script:
    input:
        trueNe = str(config["CHANGE"]["FOLDERS"]["MACRO"]) + '/forward.ne',
        map = str(config["CHANGE"]["FOLDERS"]["MACRO"]) + '/uniform.map',
    output:
        [f"{sim.FOLDER}/slimulation.slim".replace(" ","") for sim in sims.itertuples()],
    params:
        macro=str(config['CHANGE']['FOLDERS']['MACRO']),
        micro=str(config['CHANGE']['FOLDERS']['MICRO']),
        n=str(config['CHANGE']['SIMULATE']['SAMPSIZE']),
        ageSplit=str(config['CHANGE']['SIMULATE']['TSPLIT']),
        V=str(1),
        L=str(config['CHANGE']['SIMULATE']['CMLEN']),
        m=str(config['CHANGE']['SIMULATE']['NUMSUBPOP']),
        q=str(config['CHANGE']['SIMULATE']['MIGRRATE']),
        rho=str(config['CHANGE']['SIMULATE']['RHO']),
        gcProp=str(config['CHANGE']['SIMULATE']['GCPROP']),
        gcLen=str(config['CHANGE']['SIMULATE']['GCLEN']),
        a=str(config['FIXED']['SIMULATE']['a']),
        b=str(config['FIXED']['SIMULATE']['b']),
        sampleequal=str(config['FIXED']['SIMULATE']['SAMPLEEQUAL']),
    shell:
        """
        python scripts/writeSlimDemography.py \
            {params.macro} \
            {params.micro} \
            {input.trueNe} \
            {params.n} \
            {params.ageSplit} \
            {params.V} \
            {params.L} \
            {params.m} \
            {params.q} \
            {params.rho} \
            {params.gcProp} \
            {params.gcLen} \
            {params.a} \
            {params.b} \
            {params.sampleequal}
        """

### generate tree, vcf data ###

# slim forward
rule slim:
    input:
        "{macro}/{micro}/{seed}/slimulation.slim",
    output:
        trees="{macro}/{micro}/{seed}/slimulation.trees",
        freq="{macro}/{micro}/{seed}/slimulation.freq",
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}',
    shell:
        '''
        ../../software/slim \
            {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/slimulation.slim \
        '''

# msprime backward
g = open(str(config['CHANGE']['SIMULATE']['tNe']),'r')
for line in g:
    itr, ancNe = line.strip().split('\t')
rule msprime:
    input:
        trees = "{macro}/{micro}/{seed}/slimulation.trees",
    output:
        bcf="{macro}/{micro}/{seed}/slimulation.bcf.gz",
    params:
        mu=str(config['CHANGE']['SIMULATE']['MU']),
        rho=str(config['CHANGE']['SIMULATE']['RHO']),
        ancNe=str(ancNe),
        maf=str(config['CHANGE']['SIMULATE']['MSPMAF']),
        sampsize=str(config['CHANGE']['SIMULATE']['SAMPSIZE']),
        ploidy=str(ploidy),
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}',
    shell:
        """
        ancNe=$(tail -n 1 | cut -f 2)
        python scripts/treeVCF.py \
            {output.bcf} \
            {input.trees} \
            {params.mu} \
            {params.rho} \
            {params.ancNe} \
            {params.maf} \
            {params.sampsize} \
            {params.ploidy}
        """

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

# true phase
rule genotyping_error:
    input:
        vcf='{macro}/{micro}/{seed}/slimulation.vcf.gz',
    output:
        out='{macro}/{micro}/{seed}/large.chr1.vcf.gz',
    params:
        gter=str(config['CHANGE']['SIMULATE']['GTERR']),
    shell:
        """
        python ../../scripts/utilities/add-uniform-err.py \
            --input_file {input.vcf} \
            --output_file {output.out} \
            --unif_err_rate {params.gter}
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

##### haplotype phasing

### rephasing

# remove the phase in reference samples
rule unphase_ref:
    input:
        refvcf='{study}/gtdata/refpop/chr{num}.shrink.vcf.gz',
    output:
        refvcf='{study}/gtdata/refpop/chr{num}.unphased.vcf.gz',
    shell:
        '''
        python ../../scripts/utilities/remove-phase.py \
            --phased_vcf_file {input.refvcf} \
            --unphased_vcf_file {output.refvcf}.temp
        zcat {output.refvcf}.temp | bgzip -c > {output.refvcf}
        rm -f {output.refvcf}.temp
        '''

# remove phase in admixed samples
rule unphase_adx:
    input:
        adxvcf='{study}/gtdata/adxpop/chr{num}.shrink.vcf.gz',
    output:
        adxvcf='{study}/gtdata/adxpop/chr{num}.unphased.vcf.gz',
    shell:
        '''
        python ../../scripts/utilities/remove-phase.py \
            --phased_vcf_file {input.adxvcf} \
            --unphased_vcf_file {output.adxvcf}.temp
        zcat {output.adxvcf}.temp | bgzip -c > {output.adxvcf}
        rm -f {output.adxvcf}.temp
        '''

rule merge_vcfs:
    input:
        adxvcf='{study}/gtdata/adxpop/chr{num}.unphased.vcf.gz',
        refvcf='{study}/gtdata/refpop/chr{num}.unphased.vcf.gz',
    output:
        allvcf='{study}/gtdata/all/chr{num}.unphased.vcf.gz',
    shell:
        '''
        tabix -fp vcf {input.adxvcf}
        tabix -fp vcf {input.refvcf}
        bcftools query -f "%CHROM\t%POS\n" {input.adxvcf} > {input.adxvcf}.pos
        bcftools query -f "%CHROM\t%POS\n" {input.refvcf} > {input.refvcf}.pos
        mkdir -p {wildcards.study}/gtdata/all
        python ../../scripts/pre-processing/shared.py \
            {input.adxvcf}.pos \
            {input.refvcf}.pos \
            {wildcards.study}/gtdata/all/chr{wildcards.num}.intersection.Rfile.txt
        bcftools merge \
            -O z \
            -R {wildcards.study}/gtdata/all/chr{wildcards.num}.intersection.Rfile.txt \
            -o {output.allvcf} \
            {input.adxvcf} {input.refvcf}
        rm -f {input.adxvcf}.pos {input.refvcf}.pos
        '''

# another phasing strategy
# phase using ref and admixs all together
# helps w/ phase consistency
# this controls switch errors better
rule phase_all:
    input:
        allvcf='{study}/gtdata/all/chr{num}.unphased.vcf.gz',
        chrmap='{study}/maps/chr{num}.map',
    output:
        allvcf='{study}/gtdata/all/chr{num}.rephased.vcf.gz',
    params:
        allvcfout='{study}/gtdata/all/chr{num}.rephased',
        xmx=config['change']['xmxmem'],
        excludesamples=str(config['change']['existing-data']['exclude-samples']),
        window=str(config['fixed']['beagle-parameters']['window']),
        rnsd=str(config['fixed']['beagle-parameters']['random-seed']),
    shell:
        '''
        java -Xmx{params.xmx}g -jar ../../software/beagle.jar \
            gt={input.allvcf} \
            map={input.chrmap} \
            out={params.allvcfout} \
            excludesamples={params.excludesamples} \
            window={params.window} \
            seed={params.rnsd}
        '''

# subset the phased files for admixed samples
rule subset_phased_adx:
    input:
        allvcf='{study}/gtdata/all/chr{num}.rephased.vcf.gz',
        adxsam='{study}/gtdata/adxpop/chr{num}.sample.txt'
    output:
        adxvcf='{study}/gtdata/adxpop/chr{num}.rephased.vcf.gz',
    shell:
        '''
        tabix -fp vcf {input.allvcf}
        bcftools view \
            -S {input.adxsam} \
            --force-samples \
            -O z \
            -o {output.adxvcf} \
            {input.allvcf}
        rm -f {input.adxsam}
        '''

# subset the phased files for reference samples
rule subset_phased_ref:
    input:
        allvcf='{study}/gtdata/all/chr{num}.rephased.vcf.gz',
        refsam='{study}/gtdata/refpop/chr{num}.sample.txt'
    output:
        refvcf='{study}/gtdata/refpop/chr{num}.rephased.vcf.gz',
    shell:
        '''
        tabix -fp vcf {input.allvcf}
        bcftools view \
            -S {input.refsam} \
            --force-samples \
            -O z \
            -o {output.refvcf} \
            {input.allvcf}
        rm -f {input.refsam}
        '''

### initial phasing

rule phase_ref:
    input:
        adxvcf='{study}/gtdata/adxpop/chr{num}.shrink.vcf.gz',
        refvcf='{study}/gtdata/refpop/chr{num}.shrink.vcf.gz',
        chrmap='{study}/maps/chr{num}.map',
    output:
        adxvcf='{study}/gtdata/adxpop/chr{num}.referencephased.vcf.gz',
    params:
        adxvcfout='{study}/gtdata/adxpop/chr{num}.referencephased',
        xmx=config['change']['xmxmem'],
        excludesamples=str(config['change']['existing-data']['exclude-samples']),
        impute=str(config['fixed']['beagle-parameters']['impute']),
        window=str(config['fixed']['beagle-parameters']['window']),
        rnsd=str(config['fixed']['beagle-parameters']['random-seed']),
    shell:
        '''
        java -Xmx{params.xmx}g -jar ../../software/beagle.jar \
            gt={input.adxvcf} \
            ref={input.refvcf} \
            map={input.chrmap} \
            out={params.adxvcfout} \
            excludesamples={params.excludesamples} \
            impute={params.impute} \
            window={params.window} \
            seed={params.rnsd}
        '''
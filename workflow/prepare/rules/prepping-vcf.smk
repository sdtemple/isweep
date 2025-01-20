##### prepare the vcf files

### inputs, string management, count sample size, make mac
macro=str(config['change']['want-data']['your-analysis-folder'])
low=int(float(str(config['change']['want-data']['chr-low'])))
high=int(float(str(config['change']['want-data']['chr-high'])))
reffolder=str(config['change']['existing-data']['ref-folder'])
refprefix=str(config['change']['existing-data']['ref-prefix'])
refsuffix=str(config['change']['existing-data']['ref-suffix'])
adxfolder=str(config['change']['existing-data']['adx-folder'])
adxprefix=str(config['change']['existing-data']['adx-prefix'])
adxsuffix=str(config['change']['existing-data']['adx-suffix'])

### ln -s the vcf files

# direct to the vcf file for reference samples
rule copy_ln_vcf_ref:
    input:
        refvcf=reffolder + '/' + refprefix + '{num}' + refsuffix,
        maps='{study}/maps/chr{num}.map',
    output:
        refvcf='{study}/gtdata/refpop/chr{num}.vcf.gz',
    shell:
        '''
        ln -s {input.refvcf} {output.refvcf}
        '''

# direct to the vcf file for admixed samples
rule copy_ln_vcf_adx:
    input:
        adxvcf=adxfolder + '/' + adxprefix + '{num}' + adxsuffix,
        maps='{study}/maps/chr{num}.map',
    output:
        adxvcf='{study}/gtdata/adxpop/chr{num}.vcf.gz',
    shell:
        '''
        ln -s {input.adxvcf} {output.adxvcf}
        '''

### filter vcf files by minimum allele count

rule shrink_vcf_adx:
    input:
        adxvcf='{study}/gtdata/adxpop/chr{num}.vcf.gz',
    output:
        adxvcfshrink='{study}/gtdata/adxpop/chr{num}.shrink.vcf.gz',
    params:
        minmac=str(config['change']['bcftools-parameters']['c-min-mac']),
        keepsamples=str(config['change']['existing-data']['keep-samples']),
        chrnamemap=str(config['change']['existing-data']['rename-chrs-map-adx']),
    shell:
        '''
        tabix -fp vcf {input.adxvcf}
        bcftools view \
            -c {params.minmac}:nonmajor \
            -v snps \
            -S {params.keepsamples} \
            --force-samples \
            -O z \
            -o {output.adxvcfshrink}.unannotated \
            {input.adxvcf}
        tabix -fp vcf {output.adxvcfshrink}.unannotated
        bcftools annotate \
            -O z \
            -o {output.adxvcfshrink} \
            --rename-chrs {params.chrnamemap} \
            {output.adxvcfshrink}.unannotated
        rm -f {output.adxvcfshrink}.unannotated
        rm -f {output.adxvcfshrink}.unannotated.tbi
        '''

rule shrink_vcf_ref:
    input:
        refvcf='{study}/gtdata/refpop/chr{num}.vcf.gz',
    output:
        refvcfshrink='{study}/gtdata/refpop/chr{num}.shrink.vcf.gz',
    params:
        minmac=str(config['change']['bcftools-parameters']['c-min-mac']),
        chrnamemap=str(config['change']['existing-data']['rename-chrs-map-ref']),
    shell:
        '''
        tabix -fp vcf {input.refvcf}
        bcftools view \
            -c {params.minmac}:nonmajor \
            -v snps \
            -O z \
            -o {output.refvcfshrink}.unannotated \
            {input.refvcf}
        tabix -fp vcf {output.refvcfshrink}.unannotated
        bcftools annotate \
            -O z \
            -o {output.refvcfshrink} \
            --rename-chrs {params.chrnamemap} \
            {output.refvcfshrink}.unannotated
        rm -f {output.refvcfshrink}.unannotated
        rm -f {output.refvcfshrink}.unannotated.tbi
        '''

### write samples text file

# write the sample names for reference samples
rule write_ref_sample_names:
    input:
        refvcf='{study}/gtdata/refpop/chr{num}.shrink.vcf.gz',
    output:
        refsample='{study}/gtdata/refpop/chr{num}.sample.txt',
    shell:
        '''
        bcftools query -l {input.refvcf} > {output.refsample}
        '''

# write the sample names for admixed samples
rule write_adx_sample_names:
    input:
        adxvcf='{study}/gtdata/adxpop/chr{num}.shrink.vcf.gz',
    output:
        adxsample='{study}/gtdata/adxpop/chr{num}.sample.txt',
    shell:
        '''
        bcftools query -l {input.adxvcf} > {output.adxsample}
        '''
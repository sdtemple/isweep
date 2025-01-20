##### prepare the vcf files

### inputs, string management, count sample size, make mac
macro=str(config['change']['want-data']['your-analysis-folder'])
low=int(float(str(config['change']['want-data']['chr-low'])))
high=int(float(str(config['change']['want-data']['chr-high'])))
refgdsfolder=str(config['change']['existing-data']['ref-folder'])
refgdsprefix=str(config['change']['existing-data']['ref-prefix'])
refgdssuffix=str(config['change']['existing-data']['ref-suffix'])
adxgdsfolder=str(config['change']['existing-data']['adx-folder'])
adxgdsprefix=str(config['change']['existing-data']['adx-prefix'])
adxgdssuffix=str(config['change']['existing-data']['adx-suffix'])


### ln -s the gds files

# direct to the gds file for reference samples
rule copy_ln_gds_ref:
    input:
        refgds=refgdsfolder + '/' + refgdsprefix + '{num}' + refgdssuffix,
        maps='{study}/maps/chr{num}.map',
    output:
        refgds='{study}/gtdata/refpop/chr{num}.gds',
    shell:
        '''
        ln -s {input.refgds} {output.refgds}
        '''

# direct to the gds file for admixed samples
rule copy_ln_gds_adx:
    input:
        adxgds=adxgdsfolder + '/' + adxgdsprefix + '{num}' + adxgdssuffix,
        maps='{study}/maps/chr{num}.map',
    output:
        adxgds='{study}/gtdata/adxpop/chr{num}.gds',
    shell:
        '''
        ln -s {input.adxgds} {output.adxgds}
        '''

### convert gds to vcf

# convert gds to vcf for reference samples
rule gds_to_vcf_ref:
    input:
        refgds='{study}/gtdata/refpop/chr{num}.gds',
    output:
        refvcf='{study}/gtdata/refpop/chr{num}.vcf.gz',
    shell:
        '''
        Rscript --vanilla ../../scripts/pre-processing/gds-to-vcf.R {input.refgds} {output.refvcf}
        rm -f {input.refgds}.seq.gds
        '''

# convert gds to vcf for admixed samples
rule gds_to_vcf_adx:
    input:
        adxgds='{study}/gtdata/adxpop/chr{num}.gds',
    output:
        adxvcf='{study}/gtdata/adxpop/chr{num}.vcf.gz',
    params:
        keepsamples=str(config['change']['existing-data']['keep-samples']),
        minmac=str(config['change']['bcftools-parameters']['c-min-mac']),
        minmis=str(config['change']['bcftools-parameters']['missingness']),
    shell:
        '''
        Rscript --vanilla ../../scripts/pre-processing/subset-gds.R \
            {input.adxgds} \
            {params.keepsamples} \
            {params.minmac} \
            {params.minmis} \
            {wildcards.study}/gtdata/adxpop/chr{wildcards.num}
        rm -f {input.adxgds}.seq.gds
        '''

### filter vcf files by minimum allele count

rule shrink_vcf_adx:
    input:
        adxvcf='{study}/gtdata/adxpop/chr{num}.vcf.gz',
    output:
        adxvcfshrink='{study}/gtdata/adxpop/chr{num}.shrink.vcf.gz',
    params:
        keepsamples=str(config['change']['existing-data']['keep-samples']),
        minmac=str(config['change']['bcftools-parameters']['c-min-mac']),
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

rule merge_vcfs:
    input:
        adxvcf='{study}/gtdata/adxpop/chr{num}.shrink.vcf.gz',
        refvcf='{study}/gtdata/refpop/chr{num}.shrink.vcf.gz',
    output:
        allvcf='{study}/gtdata/all/chr{num}.vcf.gz',
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

# subset the phased files for admixed samples
rule subset_phased_adx:
    input:
        allvcf='{study}/gtdata/all/chr{num}.vcf.gz',
        adxsam='{study}/gtdata/adxpop/chr{num}.sample.txt'
    output:
        adxvcf='{study}/gtdata/adxpop/chr{num}.shared.vcf.gz',
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
        allvcf='{study}/gtdata/all/chr{num}.vcf.gz',
        refsam='{study}/gtdata/refpop/chr{num}.sample.txt'
    output:
        refvcf='{study}/gtdata/refpop/chr{num}.shared.vcf.gz',
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

# perform local ancestry inference
# using flare software
rule flare:
    input:
        refvcf='{study}/gtdata/refpop/chr{num}.vcf.shared.gz',
        adxvcf='{study}/gtdata/adxpop/chr{num}.vcf.shared.gz',
        allvcf='{study}/gtdata/all/chr{num}.vcf.gz',
        chrmap='{study}/maps/chr{num}.map',
        refpanelmap=str(config['change']['existing-data']['ref-panel-map']),
    output:
        outvcf='{study}/lai/chr{num}.flare.anc.vcf.gz',
    params:
        xmxmem=str(config['change']['xmxmem']),
        gen=str(config['fixed']['flare-parameters']['gen']),
        minmaf=str(config['fixed']['flare-parameters']['min-maf']),
        minmac=str(config['fixed']['flare-parameters']['min-mac']),
        probs=str(config['change']['flare-parameters']['probs']),
        out='{study}/lai/chr{num}.flare'
    shell:
        '''
        mkdir -p {wildcards.study}/lai
        java -Xmx{params.xmxmem}g -jar ../../software/flare.jar \
            ref={input.refvcf} \
            ref-panel={input.refpanelmap} \
            gt={input.adxvcf} \
            map={input.chrmap} \
            out={params.out} \
            gen={params.gen} \
            min-maf={params.minmaf} \
            min-mac={params.minmac} \
            probs={params.probs}
        rm -f {wildcards.study}/gtdata/refpop/chr{wildcards.num}.shared.vcf.gz
        rm -f {wildcards.study}/gtdata/refpop/chr{wildcards.num}.vcf.gz.tbi
        rm -f {wildcards.study}/gtdata/adxpop/chr{wildcards.num}.shared.vcf.gz
        rm -f {wildcards.study}/gtdata/adxpop/chr{wildcards.num}.vcf.gz.tbi
        rm -f {input.allvcf}
        rm -f {input.allvcf}.tbi 
        '''

# for admixed samples
rule hapibd_adx:
    input:
        adxvcf='{study}/gtdata/adxpop/chr{num}.shrink.vcf.gz',
        chrmap='{study}/maps/chr{num}.map',
    output:
        adxhap='{study}/ibdsegs/chr{num}.adx.hapibd.ibd.gz',
    params:
        minout=str(config['fixed']['hapibd-parameters']['min-output']),
        minextend=str(config['fixed']['hapibd-parameters']['min-extend']),
        minseed=str(config['fixed']['hapibd-parameters']['min-seed']),
        maxgap=str(config['fixed']['hapibd-parameters']['max-gap']),
        adxout='{study}/ibdsegs/chr{num}.adx.hapibd',
        xmxmem=str(config['change']['xmxmem']),
        minmac=str(config['fixed']['hapibd-parameters']['min-mac']),
    shell:
        '''
        mkdir -p {wildcards.study}/ibdsegs
        java -Xmx{params.xmxmem}g -jar ../../software/hap-ibd.jar \
            gt={input.adxvcf} \
            map={input.chrmap} \
            out={params.adxout} \
            min-output={params.minout} \
            min-extend={params.minextend} \
            min-seed={params.minseed} \
            max-gap={params.maxgap} \
            min-mac={params.minmac}
        '''

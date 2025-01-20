##### local ancestry inference

# perform local ancestry inference
# using flare software
rule flare_rephased:
    input:
        refvcf='{study}/gtdata/refpop/chr{num}.rephased.vcf.gz',
        adxvcf='{study}/gtdata/adxpop/chr{num}.rephased.vcf.gz',
        allvcf='{study}/gtdata/all/chr{num}.rephased.vcf.gz',
        chrmap='{study}/maps/chr{num}.map',
        refpanelmap=str(config['change']['existing-data']['ref-panel-map']),
    output:
        outvcf='{study}/lai/chr{num}.rephased.flare.anc.vcf.gz',
    params:
        software=str(config['change']['pipe']['software']),
        nthreads=str(config['change']['cluster-resources']['threads']),
        xmxmem=str(config['change']['cluster-resources']['xmxmem']),
        gen=str(config['fixed']['flare-parameters']['gen']),
        minmaf=str(config['fixed']['flare-parameters']['min-maf']),
        minmac=str(config['fixed']['flare-parameters']['min-mac']),
        probs=str(config['change']['flare-parameters']['probs']),
        prog=str(config['fixed']['programs']['flare']),
        out='{study}/lai/chr{num}.rephased.flare'
    shell:
        '''
        mkdir -p {wildcards.study}/lai
        java -Xmx{params.xmxmem}g -jar {params.software}/{params.prog} \
            ref={input.refvcf} \
            ref-panel={input.refpanelmap} \
            gt={input.adxvcf} \
            map={input.chrmap} \
            out={params.out} \
            gen={params.gen} \
            min-maf={params.minmaf} \
            min-mac={params.minmac} \
            probs={params.probs} \
            nthreads={params.nthreads}
        rm -f {wildcards.study}/gtdata/adxpop/chr{wildcards.num}.unphased.vcf.gz
        rm -f {wildcards.study}/gtdata/adxpop/chr{wildcards.num}.unphased.vcf.gz.tbi
        rm -f {wildcards.study}/gtdata/refpop/chr{wildcards.num}.unphased.vcf.gz
        rm -f {wildcards.study}/gtdata/refpop/chr{wildcards.num}.unphased.vcf.gz.tbi
        rm -f {wildcards.study}/gtdata/all/chr{wildcards.num}.unphased.vcf.gz
        rm -f {wildcards.study}/gtdata/refpop/chr{wildcards.num}.vcf.gz
        rm -f {wildcards.study}/gtdata/refpop/chr{wildcards.num}.vcf.gz.tbi
        rm -f {wildcards.study}/gtdata/adxpop/chr{wildcards.num}.vcf.gz
        rm -f {wildcards.study}/gtdata/adxpop/chr{wildcards.num}.vcf.gz.tbi
        rm -f {input.allvcf}
        rm -f {input.allvcf}.tbi 
        '''

# perform local ancestry inference
# using flare software
rule flare_reference_phased:
    input:
        refvcf='{study}/gtdata/refpop/chr{num}.shrink.vcf.gz',
        adxvcf='{study}/gtdata/adxpop/chr{num}.referencephased.vcf.gz',
        chrmap='{study}/maps/chr{num}.map',
        refpanelmap=str(config['change']['existing-data']['ref-panel-map']),
    output:
        outvcf='{study}/lai/chr{num}.referencephased.flare.anc.vcf.gz',
    params:
        software=str(config['change']['pipe']['software']),
        nthreads=str(config['change']['cluster-resources']['threads']),
        xmxmem=str(config['change']['cluster-resources']['xmxmem']),
        gen=str(config['fixed']['flare-parameters']['gen']),
        minmaf=str(config['fixed']['flare-parameters']['min-maf']),
        minmac=str(config['fixed']['flare-parameters']['min-mac']),
        probs=str(config['change']['flare-parameters']['probs']),
        prog=str(config['fixed']['programs']['flare']),
        out='{study}/lai/chr{num}.referencephased.flare'
    shell:
        '''
        mkdir -p {wildcards.study}/lai
        java -Xmx{params.xmxmem}g -jar {params.software}/{params.prog} \
            ref={input.refvcf} \
            ref-panel={input.refpanelmap} \
            gt={input.adxvcf} \
            map={input.chrmap} \
            out={params.out} \
            gen={params.gen} \
            min-maf={params.minmaf} \
            min-mac={params.minmac} \
            probs={params.probs} \
            nthreads={params.nthreads}
        rm -f {wildcards.study}/gtdata/refpop/chr{wildcards.num}.shrink.vcf.gz
        rm -f {wildcards.study}/gtdata/adxpop/chr{wildcards.num}.shrink.vcf.gz
        '''
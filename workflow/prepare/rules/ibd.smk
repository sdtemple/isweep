##### ibd segment detection

### perform ibd segment detection
### using hap-ibd software

# for admixed samples
rule hapibd_adx_referencephased:
    input:
        adxvcf='{study}/gtdata/adxpop/chr{num}.referencephased.vcf.gz',
        chrmap='{study}/maps/chr{num}.map',
    output:
        adxhap='{study}/ibdsegs/chr{num}.referencephased.adx.hapibd.ibd.gz',
    params:
        minout=str(config['fixed']['hapibd-parameters']['min-output']),
        minextend=str(config['fixed']['hapibd-parameters']['min-extend']),
        minseed=str(config['fixed']['hapibd-parameters']['min-seed']),
        maxgap=str(config['fixed']['hapibd-parameters']['max-gap']),
        adxout='{study}/ibdsegs/chr{num}.referencephased.adx.hapibd',
        xmxmem=str(config['change']['xmxmem']),
        minmac=str(config['fixed']['hapibd-parameters']['min-mac']),
    shell:
        '''
        mkdir -p {wildcards.study}/ibdsegs
        java -Xmx{params.xmxmem}g -jar ../../hap-ibd.jar \
            gt={input.adxvcf} \
            map={input.chrmap} \
            out={params.adxout} \
            min-output={params.minout} \
            min-extend={params.minextend} \
            min-seed={params.minseed} \
            max-gap={params.maxgap} \
            min-mac={params.minmac}
        '''

# for admixed samples
rule hapibd_adx_rephased:
    input:
        adxvcf='{study}/gtdata/adxpop/chr{num}.rephased.vcf.gz',
        chrmap='{study}/maps/chr{num}.map',
    output:
        adxhap='{study}/ibdsegs/chr{num}.rephased.adx.hapibd.ibd.gz',
    params:
        minout=str(config['fixed']['hapibd-parameters']['min-output']),
        minextend=str(config['fixed']['hapibd-parameters']['min-extend']),
        minseed=str(config['fixed']['hapibd-parameters']['min-seed']),
        maxgap=str(config['fixed']['hapibd-parameters']['max-gap']),
        adxout='{study}/ibdsegs/chr{num}.rephased.adx.hapibd',
        xmxmem=str(config['change']['xmxmem']),
        minmac=str(config['fixed']['hapibd-parameters']['min-mac']),
    shell:
        '''
        mkdir -p {wildcards.study}/ibdsegs
        java -Xmx{params.xmxmem}g -jar ../../hap-ibd.jar \
            gt={input.adxvcf} \
            map={input.chrmap} \
            out={params.adxout} \
            min-output={params.minout} \
            min-extend={params.minextend} \
            min-seed={params.minseed} \
            max-gap={params.maxgap} \
            min-mac={params.minmac}
        '''

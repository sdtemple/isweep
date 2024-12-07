# make excludesamples file for browning lab jar

# inputs, string management
vcffolder=str(config['CHANGE']['EXISTING']['VCFS'])
vcfpre=str(config['CHANGE']['EXISTING']['VCFPRE'])
vcfsuf=str(config['CHANGE']['EXISTING']['VCFSUF'])
low=str(config['CHANGE']['ISWEEP']['CHRHIGH'])
high=str(config['CHANGE']['ISWEEP']['CHRLOW'])
subsample=str(config['CHANGE']['ISWEEP']['SUBSAMPLE'])

# get all samples from bcftools
rule make_samples:
    input:
        vcf=vcffolder + '/' + vcfpre + high + vcfsuf,
    output:
        sample=macro+'/sample.txt',
    shell:
        'bcftools query -l {input.vcf} > {output.sample}'

# copy input subsample to subsample.txt
rule make_subsample:
    input:
        filein=str(subsample),
    output:
        fileout=macro+'/subsample.txt',
    shell:
        'cp {input.filein} {output.fileout}'

# make excludesamples
rule make_excludesamples:
    input:
        sample=macro+'/sample.txt',
        subsample=macro+'/subsample.txt',
    output:
        exclsample=macro+'/excludesamples.txt',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
    shell:
        """
        python {params.scripts}/exclude-samples.py \
            {input.sample} \
            {input.subsample} \
            {output.exclsample}
        """

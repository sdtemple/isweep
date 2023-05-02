# make excludesamples file for Browning JAR files

vcffolder=str(config['CHANGE']['EXISTING']['VCFS'])
vcfpre=str(config['CHANGE']['EXISTING']['VCFPRE'])
vcfsuf=str(config['CHANGE']['EXISTING']['VCFSUF'])
low=str(config['CHANGE']['ISWEEP']['CHRHIGH'])
high=str(config['CHANGE']['ISWEEP']['CHRLOW'])
subsample=str(config['CHANGE']['ISWEEP']['SUBSAMPLE'])

rule make_samples:
    input:
        vcf=vcffolder + '/' + vcfpre + high + vcfsuf,
    output:
        sample=macro+'/sample.txt',
    shell:
        'bcftools query -l {input.vcf} > {output.sample}'

rule make_subsample:
    input:
        filein=str(subsample),
    output:
        fileout=macro+'/subsample.txt',
    shell:
        'cp {input.filein} {output.fileout}'

rule make_excludesamples:
    input:
        sample=macro+'/sample.txt',
        subsample=macro+'/subsample.txt',
    output:
        exclsample=macro+'/excludesamples.txt',
    shell:
        'python {config[CHANGE][FOLDERS][TERMINALSCRIPTS]}/exclude-samples.py {input.sample} {input.subsample} {output.exclsample}'

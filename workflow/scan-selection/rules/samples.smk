# Reformat the sample files, which includes
# flle for excludesamples= in Browning Lab software.

localrules: make_samples, make_subsample

# get all samples from bcftools
rule make_samples:
    input:
        vcf=vcffolder + '/' + vcfpre + str(high) + vcfsuf,
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
    shell:
        """
        python ../../scripts/utilities/exclude-samples.py \
            --sample_file {input.sample} \
            --subsample_file {input.subsample} \
            --exclude_sample_file {output.exclsample}
        """

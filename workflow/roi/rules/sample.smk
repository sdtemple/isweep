# make excludesamples file for Browning JAR files

rule make_samples:
    input:
        vcf='{config[CHANGE][EXISTING][VCFS]}/{config[CHANGE][EXISTING][VCFPRE]}{config[CHANGE][ISWEEP][CHRLOW]}{config[CHANGE][EXISTING][VCFSUF]}',
    output:
        sample='{config[CHANGE][FOLDERS][STUDY]}/sample.txt',
    shell:
        'bcftools query -l {input.vcf} > {output.sample}'

rule make_excludesamples:
    input:
        samples='{config[CHANGE][FOLDERS][STUDY]}/sample.txt',
        subsample='{config[CHANGE][ISWEEP][SUBSAMPLE]}',
    output:
        '{config[CHANGE][FOLDERS][STUDY]}/excludesamples.txt',
    script:
        '{config[CHANGE][FOLDERS][SNAKESCRIPTS]}/exclude-samples.py'

rule make_subsample:
    input:
        samples='{config[CHANGE][FOLDERS][STUDY]}/sample.txt',
        subsample='{config[CHANGE][FOLDERS][STUDY]}/excludesamples.txt',
    output:
        '{config[CHANGE][FOLDERS][STUDY]}/subsample.txt',
    script:
        '{config[CHANGE][FOLDERS][SNAKESCRIPTS]}/exclude-samples.py'

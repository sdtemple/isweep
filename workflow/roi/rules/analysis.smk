# conduct analysis using isweep package
# seth temple, sdtemple@uw.edu
# april 29, 2023

rule rank:
    input:
        short='{cohort}/roi{roi}/chr{chr}/center{center}/left{left}/right{right}/focus.ibd.gz',
        vcf='{cohort}/roi{roi}/chr{chr}/center{center}/left{left}/right{right}/chr{chr}.focused.vcf.gz',
    output:
        fileout='{cohort}/roi{roi}/chr{chr}/center{center}/left{left}/right{right}/isweep.ranks.tsv.gz',
    params:
        diameter='{config[FIXED][ISWEEP][DIAMETER]}',
        atleast='{config[FIXED][ISWEEP][ATLEAST]}',
        q1='{config[FIXED][ISWEEP][MINAF]}',
        q2='{config[FIXED][ISWEEP][MAXAF]}',
        rulesigma='{config[FIXED][ISWEEP][RULESIGMA]}',
    shell:
        'python {config[CHANGE][FOLDERS][TERMINALSCRIPTS]}/rank-isweep.py {input.short} {input.vcf} {output.fileout} {params.diameter} {params.atleast} {params.q1} {params.q2} {params.rulesigma}'

rule infer:
    input:
        long='{cohort}/{roi}/chr{chr}/center{center}/left{left}/right{right}/mom.ibd.gz',
        fileout='{cohort}/roi{roi}/chr{chr}/center{center}/left{left}/right{right}/isweep.ranks.tsv.gz',
    output:
        fileout='{cohort}/roi{roi}/chr{chr}/center{center}/left{left}/right{right}/isweep.inferences.tsv.gz',
    params:
        nboot='{config[FIXED][ISWEEP][NBOOT]}',
        cm='{config[FIXED][ISWEEP][MOMCUTOFF]}',
        n='{config[CHANGE][ISWEEP][SAMPSIZE]}',
        effdemo='{config[CHANGE][ISWEEP][NE]}',
        ploidy='{config[FIXED][HAPIBD][PLOIDY]}',
    shell:
        'python {config[CHANGE][FOLDERS][TERMINALSCRIPTS]}/infer-isweep.py {input.short} {input.vcf} {output.fileout} {params.diameter} {params.atleast} {params.q1} {params.q2} {params.rulesigma}; touch {config[CHANGE][FOLDERS][STUDY]}/done{roi}'

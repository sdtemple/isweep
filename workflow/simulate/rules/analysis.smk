# conduct analysis using isweep package
# seth temple, sdtemple@uw.edu
# may 3, 2023

# rank polymorphisms in focus region
rule rank:
    input:
        short='{macro}/{micro}/{seed}/short.chr1.ibd.gz',
        vcf='{macro}/{micro}/{seed}/small.chr1.vcf.gz',
    output:
        fileout='{macro}/{micro}/{seed}/isweep.ranks.tsv.gz',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        diameter=str(config['FIXED']['ISWEEP']['DIAMETER']),
        atleast=str(config['FIXED']['ISWEEP']['ATLEAST']),
        q1=str(config['FIXED']['ISWEEP']['MINAAF']),
        q2=str(config['FIXED']['ISWEEP']['MAXAAF']),
        rulesigma=str(config['FIXED']['ISWEEP']['RULESIGMA']),
    shell:
        """
        python {params.scripts}/rank-isweep.py \
            {input.short} \
            {input.vcf} \
            {output.fileout} \
            {params.diameter} \
            {params.atleast} \
            {params.q1} \
            {params.q2} \
            {params.rulesigma}
        """

rule infer:
    input:
        long='{macro}/{micro}/{seed}/long.chr1.ibd.gz',
        ranks='{macro}/{micro}/{seed}/isweep.ranks.tsv.gz',
    output:
        fileout='{macro}/{micro}/{seed}/isweep.inference.tsv.gz',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        nboot=str(config['FIXED']['ISWEEP']['NBOOT']),
        cm=str(config['FIXED']['ISWEEP']['MOMCUTOFF']),
        n=str(config['CHANGE']['SIMULATE']['SAMPSIZE']),
        ploidy=str(config['FIXED']['HAPIBD']['PLOIDY']),
        effdemo=str(config['CHANGE']['SIMULATE']['iNe']),
    shell:
        """
        python {params.scripts}/infer-isweep.py \
            {input.long} \
            {input.ranks} \
            {output.fileout} \
            {params.nboot} \
            {params.cm} \
            {params.n} \
            {params.effdemo} \
            {params.ploidy}
        """

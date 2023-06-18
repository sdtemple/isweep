# conduct analysis using isweep package
# seth temple, sdtemple@uw.edu
# may 3, 2023
# may 8, 2023 (updated with relate, clues)

# less basepairs in ranking
rule rank_vcf:
    input:
        vcfin='{macro}/{micro}/{seed}/small.chr1.vcf.gz',
    output:
        vcfout='{macro}/{micro}/{seed}/short.chr1.vcf.gz',
    params:
        maf=str(config['FIXED']['ISWEEP']['MINMAXAAF']),
    shell:
        """
        tabix -fp vcf {input.vcfin}
        bcftools view -q {params.maf}:nonmajor \
            -Oz -o {output.vcfout} \
            {input.vcfin}
        tabix -fp vcf {output.vcfout}
        """

# rank polymorphisms in focus region
rule rank:
    input:
        short='{macro}/{micro}/{seed}/short.chr1.ibd.gz',
        vcf='{macro}/{micro}/{seed}/short.chr1.vcf.gz',
    output:
        fileout='{macro}/{micro}/{seed}/isweep.ranks.tsv.gz',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        diameter=str(config['FIXED']['ISWEEP']['DIAMETER']),
        q1=str(config['FIXED']['ISWEEP']['MINMAXAAF']),
        rulesigma=str(config['FIXED']['ISWEEP']['RULESIGMA']),
    shell:
        """
        python {params.scripts}/rank-isweep.py \
            {input.short} \
            {input.vcf} \
            {output.fileout} \
            {params.diameter} \
            {params.q1} \
            {params.rulesigma}
        """

# check rank of the true polymorphism
rule rank_true:
    input:
        filein='{macro}/{micro}/{seed}/isweep.ranks.tsv.gz',
    output:
        fileout='{macro}/{micro}/{seed}/isweep.rank.true.tsv',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        loc=str(config['FIXED']['SIMULATE']['LOC']),
    shell:
        """
        python {params.scripts}/rank-true.py \
            {input.filein} \
            {output.fileout} \
            {params.loc} \
            1
        """

rule infer:
    input:
        long='{macro}/{micro}/{seed}/long.chr1.ibd.gz',
        ranks='{macro}/{micro}/{seed}/isweep.ranks.tsv.gz',
    output:
        fileout='{macro}/{micro}/{seed}/isweep.inference.tsv',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        nboot=str(config['FIXED']['ISWEEP']['NBOOT']),
        cm=str(config['FIXED']['ISWEEP']['MOMCUTOFF']),
        n=str(config['CHANGE']['SIMULATE']['SAMPSIZE']),
        ploidy=str(config['FIXED']['HAPIBD']['PLOIDY']),
        effdemo=str(config['CHANGE']['SIMULATE']['iNe']),
        quant=str(config['FIXED']['ISWEEP']['QUANT']),
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
            {params.ploidy} \
            {params.quant}
        """

##### haplotype analysis #####

rule haplotypes:
    input:
        rankin='{macro}/{micro}/{seed}/isweep.ranks.mean.tsv.gz',
        ibdwin='{macro}/{micro}/{seed}/scan.chr1.windowed.tsv.gz',
    output:
        lociout='{macro}/{micro}/{seed}/isweep.hap.pos.txt',
        freqout='{macro}/{micro}/{seed}/isweep.hap.freq.txt',
        happng='{macro}/{micro}/{seed}/isweep.hap.png',
        snppng='{macro}/{micro}/{seed}/isweep.snp.png',
    params:
        windowsize=str(config['FIXED']['ISWEEP']['WINSIZE']),
        windowstep=str(config['FIXED']['ISWEEP']['WINSTEP']),
        windowcol=str(config['FIXED']['ISWEEP']['WINCOL']),
        freqsize=str(config['FIXED']['ISWEEP']['FREQSIZE']),
        freqstep=str(config['FIXED']['ISWEEP']['FREQSTEP']),
        freqcol=str(config['FIXED']['ISWEEP']['FREQCOL']),
        scorecol=str(config['FIXED']['ISWEEP']['SCORECOL']),
        numsnp=str(config['FIXED']['ISWEEP']['NUMSNP']),
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
    shell:
        """
        python {params.scripts}/haplotypes.py \
            {input.rankin} \
            {input.ibdwin} \
            {wildcards.macro}/{wildcards.micro}/{wildcards.seed} \
            0 \
            BPWINDOW \
            CMWINDOW \
            {params.freqsize} \
            {params.freqstep} \
            {params.freqcol} \
            {params.windowsize} \
            {params.windowstep} \
            {params.windowcol} \
            {params.scorecol} \
            {params.numsnp}
        """

rule ibd_hap_pos:
    input:
        ibd='{macro}/{micro}/{seed}/scan.chr1.ibd.gz',
        locus='{macro}/{micro}/{seed}/isweep.hap.pos.txt',
    output:
        ibd='{macro}/{micro}/{seed}/long.chr1.hap.pos.ibd.gz',
    params:
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['FILTER']),
        script=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS'])+'/first-line.py',
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}'
    shell:
        """
        themean=$(python {params.script} {input.locus})
        zcat {input.ibd} | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 6 0.00 ${{themean}} | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 7 ${{themean}} 10000000000 | \
            gzip > {output.ibd}
        """

rule infer_hap:
    input:
        long='{macro}/{micro}/{seed}/long.chr1.hap.pos.ibd.gz',
        freq='{macro}/{micro}/{seed}/isweep.hap.freq.txt',
    output:
        fileout='{macro}/{micro}/{seed}/isweep.hap.tsv',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        nboot=str(config['FIXED']['ISWEEP']['NBOOT']),
        cm=str(config['FIXED']['ISWEEP']['MOMCUTOFF']),
        n=str(config['CHANGE']['SIMULATE']['SAMPSIZE']),
        ploidy=str(config['FIXED']['HAPIBD']['PLOIDY']),
        effdemo=str(config['CHANGE']['SIMULATE']['iNe']),
    shell:
        """
        freqest=$(python {params.scripts}/first-line.py {input.freq})
        python {params.scripts}/infer-hap.py \
            {input.long} \
            {output.fileout} \
            ${{freqest}} \
            {params.nboot} \
            {params.cm} \
            {params.n} \
            {params.effdemo} \
            {params.ploidy}
        """

##### refining the locus #####
# 
# rule locus:
#     input:
#         rankin='{macro}/{micro}/{seed}/isweep.ranks.mean.tsv.gz',
#         ibdwin='{macro}/{micro}/{seed}/scan.chr1.windowed.tsv.gz',
#     output:
#         lociout='{macro}/{micro}/{seed}/isweep.locus.txt',
#         freqout='{macro}/{micro}/{seed}/isweep.freq.txt',
#     params:
#         windowsize=str(config['FIXED']['ISWEEP']['WINSIZE']),
#         windowstep=str(config['FIXED']['ISWEEP']['WINSTEP']),
#         lowq=str(config['FIXED']['ISWEEP']['LOWQ']),
#         lowp=str(config['FIXED']['ISWEEP']['LOWP']),
#         scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
#     shell:
#         """
#         python {params.scripts}/locus.py \
#             {input.rankin} \
#             {input.ibdwin} \
#             {wildcards.macro}/{wildcards.micro}/{wildcards.seed} \
#             {params.windowsize} \
#             {params.windowstep} \
#             {params.lowq} \
#             {params.lowp}
#         """
#
# rule ibd_locus:
#     input:
#         ibd='{macro}/{micro}/{seed}/scan.chr1.ibd.gz',
#         locus='{macro}/{micro}/{seed}/isweep.locus.txt',
#     output:
#         ibd='{macro}/{micro}/{seed}/long.chr1.locus.ibd.gz',
#     params:
#         soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
#         prog=str(config['CHANGE']['PROGRAMS']['FILTER']),
#         script=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS'])+'/first-line.py',
#     resources:
#         mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}'
#     shell:
#         """
#         themean=$(python {params.script} {input.locus})
#         zcat {input.ibd} | \
#             java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
#             "I" 6 0.00 ${{themean}} | \
#             java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
#             "I" 7 ${{themean}} 10000000000 | \
#             gzip > {output.ibd}
#         """
#
# rule infer_locus:
#     input:
#         long='{macro}/{micro}/{seed}/long.chr1.locus.ibd.gz',
#         freq='{macro}/{micro}/{seed}/isweep.freq.txt',
#     output:
#         fileout='{macro}/{micro}/{seed}/isweep.locus.inference.tsv',
#     params:
#         scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
#         nboot=str(config['FIXED']['ISWEEP']['NBOOT']),
#         cm=str(config['FIXED']['ISWEEP']['MOMCUTOFF']),
#         n=str(config['CHANGE']['SIMULATE']['SAMPSIZE']),
#         ploidy=str(config['FIXED']['HAPIBD']['PLOIDY']),
#         effdemo=str(config['CHANGE']['SIMULATE']['iNe']),
#     shell:
#         """
#         freqest=$(python {params.scripts}/first-line.py {input.freq})
#         python {params.scripts}/infer-locus.py \
#             {input.long} \
#             {output.fileout} \
#             ${{freqest}} \
#             {params.nboot} \
#             {params.cm} \
#             {params.n} \
#             {params.effdemo} \
#             {params.ploidy}
#         """

##### decrease in coefficient of variation #####

##### probably delete this #####

rule coef_var:
    input:
        rankin='{macro}/{micro}/{seed}/isweep.ranks.mean.tsv.gz',
    output:
        coefout='{macro}/{micro}/{seed}/isweep.coef.var.txt',
    params:
        windowsize=str(config['FIXED']['ISWEEP']['WINSIZE']),
        windowstep=str(config['FIXED']['ISWEEP']['WINSTEP']),
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
    shell:
        """
        python {params.scripts}/coef-var-isweep.py \
            {input.rankin} \
            {wildcards.macro}/{wildcards.micro}/{wildcards.seed} \
            {params.windowsize} \
            {params.windowstep}
        """

rule long_ibd_coef_var:
    input:
        ibd='{macro}/{micro}/{seed}/scan.chr1.ibd.gz',
        ibdwin='{macro}/{micro}/{seed}/scan.chr1.windowed.tsv.gz',
        coefin='{macro}/{micro}/{seed}/isweep.coef.var.txt',
    output:
        ibd='{macro}/{micro}/{seed}/long.chr1.coef.var.ibd.gz',
    params:
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['FILTER']),
        script=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS'])+'/grab-freq.py',
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}'
    shell:
        """
        themean=$(python {params.script} {input.coefin})
        zcat {input.ibd} | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 6 0.00 ${{themean}} | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 7 ${{themean}} 10000000000 | \
            gzip > {output.ibd}
        """

rule infer_coef_var:
    input:
        long='{macro}/{micro}/{seed}/long.chr1.coef.var.ibd.gz',
        ranks='{macro}/{micro}/{seed}/isweep.coef.var.txt',
    output:
        fileout='{macro}/{micro}/{seed}/isweep.inference.coef.var.tsv',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        nboot=str(config['FIXED']['ISWEEP']['NBOOT']),
        cm=str(config['FIXED']['ISWEEP']['MOMCUTOFF']),
        n=str(config['CHANGE']['SIMULATE']['SAMPSIZE']),
        ploidy=str(config['FIXED']['HAPIBD']['PLOIDY']),
        effdemo=str(config['CHANGE']['SIMULATE']['iNe']),
    shell:
        """
        python {params.scripts}/infer-coef-var-isweep.py \
            {input.long} \
            {input.ranks} \
            {output.fileout} \
            {params.nboot} \
            {params.cm} \
            {params.n} \
            {params.effdemo} \
            {params.ploidy}
        """

##### everything below this is for mean, median, mode versus true loc #####

##### probably delete much of this #####

# less basepairs in ranking
rule rank_vcf_mean:
    input:
        vcfin='{macro}/{micro}/{seed}/small.chr1.mean.vcf.gz',
    output:
        vcfout='{macro}/{micro}/{seed}/short.chr1.mean.vcf.gz',
    params:
        maf=str(config['FIXED']['ISWEEP']['MINMAXAAF']),
        folder='{macro}/{micro}/{seed}',
    shell:
        """
        gunzip -c {input.vcfin} | bgzip  > {params.folder}/chrtemp.mean.vcf.bgz
        tabix -fp vcf {params.folder}/chrtemp.mean.vcf.bgz
        bcftools view -q {params.maf}:nonmajor \
            -Oz -o {output.vcfout} \
            {params.folder}/chrtemp.mean.vcf.bgz
        tabix -fp vcf {output.vcfout}
        rm {params.folder}/chrtemp.mean.vcf.bgz
        """

# rank polymorphisms in focus region
rule rank_mean:
    input:
        short='{macro}/{micro}/{seed}/short.chr1.mean.ibd.gz',
        vcf='{macro}/{micro}/{seed}/short.chr1.mean.vcf.gz',
    output:
        fileout='{macro}/{micro}/{seed}/isweep.ranks.mean.tsv.gz',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        diameter=str(config['FIXED']['ISWEEP']['DIAMETER']),
        q1=str(config['FIXED']['ISWEEP']['MINMAXAAF']),
        rulesigma=str(config['FIXED']['ISWEEP']['RULESIGMA']),
    shell:
        """
        python {params.scripts}/rank-isweep.py \
            {input.short} \
            {input.vcf} \
            {output.fileout} \
            {params.diameter} \
            {params.q1} \
            {params.rulesigma}
        """

# check rank of the true polymorphism
rule rank_true_mean:
    input:
        filein='{macro}/{micro}/{seed}/isweep.ranks.mean.tsv.gz',
    output:
        fileout='{macro}/{micro}/{seed}/isweep.rank.true.mean.tsv',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        loc=str(config['FIXED']['SIMULATE']['LOC']),
    shell:
        """
        python {params.scripts}/rank-true.py \
            {input.filein} \
            {output.fileout} \
            {params.loc} \
            1
        """

# less basepairs in ranking
rule rank_vcf_mode:
    input:
        vcfin='{macro}/{micro}/{seed}/small.chr1.mode.vcf.gz',
    output:
        vcfout='{macro}/{micro}/{seed}/short.chr1.mode.vcf.gz',
    params:
        maf=str(config['FIXED']['ISWEEP']['MINMAXAAF']),
        folder='{macro}/{micro}/{seed}',
    shell:
        """
        gunzip -c {input.vcfin} | bgzip  > {params.folder}/chrtemp.mode.vcf.bgz
        tabix -fp vcf {params.folder}/chrtemp.mode.vcf.bgz
        bcftools view -q {params.maf}:nonmajor \
            -Oz -o {output.vcfout} \
            {params.folder}/chrtemp.mode.vcf.bgz
        tabix -fp vcf {output.vcfout}
        rm {params.folder}/chrtemp.mode.vcf.bgz
        """

# rank polymorphisms in focus region
rule rank_mode:
    input:
        short='{macro}/{micro}/{seed}/short.chr1.mode.ibd.gz',
        vcf='{macro}/{micro}/{seed}/short.chr1.mode.vcf.gz',
    output:
        fileout='{macro}/{micro}/{seed}/isweep.ranks.mode.tsv.gz',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        diameter=str(config['FIXED']['ISWEEP']['DIAMETER']),
        q1=str(config['FIXED']['ISWEEP']['MINMAXAAF']),
        rulesigma=str(config['FIXED']['ISWEEP']['RULESIGMA']),
    shell:
        """
        python {params.scripts}/rank-isweep.py \
            {input.short} \
            {input.vcf} \
            {output.fileout} \
            {params.diameter} \
            {params.q1} \
            {params.rulesigma}
        """

# check rank of the true polymorphism
rule rank_true_mode:
    input:
        filein='{macro}/{micro}/{seed}/isweep.ranks.mode.tsv.gz',
    output:
        fileout='{macro}/{micro}/{seed}/isweep.rank.true.mode.tsv',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        loc=str(config['FIXED']['SIMULATE']['LOC']),
    shell:
        """
        python {params.scripts}/rank-true.py \
            {input.filein} \
            {output.fileout} \
            {params.loc} \
            1
        """

# less basepairs in ranking
rule rank_vcf_median:
    input:
        vcfin='{macro}/{micro}/{seed}/small.chr1.median.vcf.gz',
    output:
        vcfout='{macro}/{micro}/{seed}/short.chr1.median.vcf.gz',
    params:
        maf=str(config['FIXED']['ISWEEP']['MINMAXAAF']),
        folder='{macro}/{micro}/{seed}',
    shell:
        """
        gunzip -c {input.vcfin} | bgzip  > {params.folder}/chrtemp.median.vcf.bgz
        tabix -fp vcf {params.folder}/chrtemp.median.vcf.bgz
        bcftools view -q {params.maf}:nonmajor \
            -Oz -o {output.vcfout} \
            {params.folder}/chrtemp.median.vcf.bgz
        tabix -fp vcf {output.vcfout}
        rm {params.folder}/chrtemp.median.vcf.bgz
        """

# rank polymorphisms in focus region
rule rank_median:
    input:
        short='{macro}/{micro}/{seed}/short.chr1.median.ibd.gz',
        vcf='{macro}/{micro}/{seed}/short.chr1.median.vcf.gz',
    output:
        fileout='{macro}/{micro}/{seed}/isweep.ranks.median.tsv.gz',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        diameter=str(config['FIXED']['ISWEEP']['DIAMETER']),
        q1=str(config['FIXED']['ISWEEP']['MINMAXAAF']),
        rulesigma=str(config['FIXED']['ISWEEP']['RULESIGMA']),
    shell:
        """
        python {params.scripts}/rank-isweep.py \
            {input.short} \
            {input.vcf} \
            {output.fileout} \
            {params.diameter} \
            {params.q1} \
            {params.rulesigma}
        """

# check rank of the true polymorphism
rule rank_true_median:
    input:
        filein='{macro}/{micro}/{seed}/isweep.ranks.median.tsv.gz',
    output:
        fileout='{macro}/{micro}/{seed}/isweep.rank.true.median.tsv',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        loc=str(config['FIXED']['SIMULATE']['LOC']),
    shell:
        """
        python {params.scripts}/rank-true.py \
            {input.filein} \
            {output.fileout} \
            {params.loc} \
            1
        """

rule infer_mean:
    input:
        long='{macro}/{micro}/{seed}/long.chr1.mean.ibd.gz',
        ranks='{macro}/{micro}/{seed}/isweep.ranks.mean.tsv.gz',
    output:
        fileout='{macro}/{micro}/{seed}/isweep.inference.mean.tsv',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        nboot=str(config['FIXED']['ISWEEP']['NBOOT']),
        cm=str(config['FIXED']['ISWEEP']['MOMCUTOFF']),
        n=str(config['CHANGE']['SIMULATE']['SAMPSIZE']),
        ploidy=str(config['FIXED']['HAPIBD']['PLOIDY']),
        effdemo=str(config['CHANGE']['SIMULATE']['iNe']),
        quant=str(config['FIXED']['ISWEEP']['QUANT']),
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
            {params.ploidy} \
            {params.quant}
        """

rule infer_median:
    input:
        long='{macro}/{micro}/{seed}/long.chr1.median.ibd.gz',
        ranks='{macro}/{micro}/{seed}/isweep.ranks.median.tsv.gz',
    output:
        fileout='{macro}/{micro}/{seed}/isweep.inference.median.tsv',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        nboot=str(config['FIXED']['ISWEEP']['NBOOT']),
        cm=str(config['FIXED']['ISWEEP']['MOMCUTOFF']),
        n=str(config['CHANGE']['SIMULATE']['SAMPSIZE']),
        ploidy=str(config['FIXED']['HAPIBD']['PLOIDY']),
        effdemo=str(config['CHANGE']['SIMULATE']['iNe']),
        quant=str(config['FIXED']['ISWEEP']['QUANT']),
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
            {params.ploidy} \
            {params.quant}
        """

rule infer_mode:
    input:
        long='{macro}/{micro}/{seed}/long.chr1.mode.ibd.gz',
        ranks='{macro}/{micro}/{seed}/isweep.ranks.mode.tsv.gz',
    output:
        fileout='{macro}/{micro}/{seed}/isweep.inference.mode.tsv',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        nboot=str(config['FIXED']['ISWEEP']['NBOOT']),
        cm=str(config['FIXED']['ISWEEP']['MOMCUTOFF']),
        n=str(config['CHANGE']['SIMULATE']['SAMPSIZE']),
        ploidy=str(config['FIXED']['HAPIBD']['PLOIDY']),
        effdemo=str(config['CHANGE']['SIMULATE']['iNe']),
        quant=str(config['FIXED']['ISWEEP']['QUANT']),
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
            {params.ploidy} \
            {params.quant}
        """

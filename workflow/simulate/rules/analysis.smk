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

# ##### everything below this is for mean, median, mode versus true loc #####
#
# # less basepairs in ranking
# rule rank_vcf_mean:
#     input:
#         vcfin='{macro}/{micro}/{seed}/small.chr1.mean.vcf.gz',
#     output:
#         vcfout='{macro}/{micro}/{seed}/short.chr1.mean.vcf.gz',
#     params:
#         maf=str(config['FIXED']['ISWEEP']['MINMAXAAF']),
#         folder='{macro}/{micro}/{seed}',
#     shell:
#         """
#         gunzip -c {input.vcfin} | bgzip  > {params.folder}/chrtemp.mean.vcf.bgz
#         tabix -fp vcf {params.folder}/chrtemp.mean.vcf.bgz
#         bcftools view -q {params.maf}:nonmajor \
#             -Oz -o {output.vcfout} \
#             {params.folder}/chrtemp.mean.vcf.bgz
#         tabix -fp vcf {output.vcfout}
#         rm {params.folder}/chrtemp.mean.vcf.bgz
#         """
#
# # rank polymorphisms in focus region
# rule rank_mean:
#     input:
#         short='{macro}/{micro}/{seed}/short.chr1.mean.ibd.gz',
#         vcf='{macro}/{micro}/{seed}/short.chr1.mean.vcf.gz',
#     output:
#         fileout='{macro}/{micro}/{seed}/isweep.ranks.mean.tsv.gz',
#     params:
#         scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
#         diameter=str(config['FIXED']['ISWEEP']['DIAMETER']),
#         q1=str(config['FIXED']['ISWEEP']['MINMAXAAF']),
#         rulesigma=str(config['FIXED']['ISWEEP']['RULESIGMA']),
#     shell:
#         """
#         python {params.scripts}/rank-isweep.py \
#             {input.short} \
#             {input.vcf} \
#             {output.fileout} \
#             {params.diameter} \
#             {params.q1} \
#             {params.rulesigma}
#         """
#
# # check rank of the true polymorphism
# rule rank_true_mean:
#     input:
#         filein='{macro}/{micro}/{seed}/isweep.ranks.mean.tsv.gz',
#     output:
#         fileout='{macro}/{micro}/{seed}/isweep.rank.true.mean.tsv',
#     params:
#         scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
#         loc=str(config['FIXED']['SIMULATE']['LOC']),
#     shell:
#         """
#         python {params.scripts}/rank-true.py \
#             {input.filein} \
#             {output.fileout} \
#             {params.loc} \
#             1
#         """
#
# # less basepairs in ranking
# rule rank_vcf_mode:
#     input:
#         vcfin='{macro}/{micro}/{seed}/small.chr1.mode.vcf.gz',
#     output:
#         vcfout='{macro}/{micro}/{seed}/short.chr1.mode.vcf.gz',
#     params:
#         maf=str(config['FIXED']['ISWEEP']['MINMAXAAF']),
#         folder='{macro}/{micro}/{seed}',
#     shell:
#         """
#         gunzip -c {input.vcfin} | bgzip  > {params.folder}/chrtemp.mode.vcf.bgz
#         tabix -fp vcf {params.folder}/chrtemp.mode.vcf.bgz
#         bcftools view -q {params.maf}:nonmajor \
#             -Oz -o {output.vcfout} \
#             {params.folder}/chrtemp.mode.vcf.bgz
#         tabix -fp vcf {output.vcfout}
#         rm {params.folder}/chrtemp.mode.vcf.bgz
#         """
#
# # rank polymorphisms in focus region
# rule rank_mode:
#     input:
#         short='{macro}/{micro}/{seed}/short.chr1.mode.ibd.gz',
#         vcf='{macro}/{micro}/{seed}/short.chr1.mode.vcf.gz',
#     output:
#         fileout='{macro}/{micro}/{seed}/isweep.ranks.mode.tsv.gz',
#     params:
#         scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
#         diameter=str(config['FIXED']['ISWEEP']['DIAMETER']),
#         q1=str(config['FIXED']['ISWEEP']['MINMAXAAF']),
#         rulesigma=str(config['FIXED']['ISWEEP']['RULESIGMA']),
#     shell:
#         """
#         python {params.scripts}/rank-isweep.py \
#             {input.short} \
#             {input.vcf} \
#             {output.fileout} \
#             {params.diameter} \
#             {params.q1} \
#             {params.rulesigma}
#         """
#
# # check rank of the true polymorphism
# rule rank_true_mode:
#     input:
#         filein='{macro}/{micro}/{seed}/isweep.ranks.mode.tsv.gz',
#     output:
#         fileout='{macro}/{micro}/{seed}/isweep.rank.true.mode.tsv',
#     params:
#         scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
#         loc=str(config['FIXED']['SIMULATE']['LOC']),
#     shell:
#         """
#         python {params.scripts}/rank-true.py \
#             {input.filein} \
#             {output.fileout} \
#             {params.loc} \
#             1
#         """
#
# # less basepairs in ranking
# rule rank_vcf_median:
#     input:
#         vcfin='{macro}/{micro}/{seed}/small.chr1.median.vcf.gz',
#     output:
#         vcfout='{macro}/{micro}/{seed}/short.chr1.median.vcf.gz',
#     params:
#         maf=str(config['FIXED']['ISWEEP']['MINMAXAAF']),
#         folder='{macro}/{micro}/{seed}',
#     shell:
#         """
#         gunzip -c {input.vcfin} | bgzip  > {params.folder}/chrtemp.median.vcf.bgz
#         tabix -fp vcf {params.folder}/chrtemp.median.vcf.bgz
#         bcftools view -q {params.maf}:nonmajor \
#             -Oz -o {output.vcfout} \
#             {params.folder}/chrtemp.median.vcf.bgz
#         tabix -fp vcf {output.vcfout}
#         rm {params.folder}/chrtemp.median.vcf.bgz
#         """
#
# # rank polymorphisms in focus region
# rule rank_median:
#     input:
#         short='{macro}/{micro}/{seed}/short.chr1.median.ibd.gz',
#         vcf='{macro}/{micro}/{seed}/short.chr1.median.vcf.gz',
#     output:
#         fileout='{macro}/{micro}/{seed}/isweep.ranks.median.tsv.gz',
#     params:
#         scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
#         diameter=str(config['FIXED']['ISWEEP']['DIAMETER']),
#         q1=str(config['FIXED']['ISWEEP']['MINMAXAAF']),
#         rulesigma=str(config['FIXED']['ISWEEP']['RULESIGMA']),
#     shell:
#         """
#         python {params.scripts}/rank-isweep.py \
#             {input.short} \
#             {input.vcf} \
#             {output.fileout} \
#             {params.diameter} \
#             {params.q1} \
#             {params.rulesigma}
#         """
#
# # check rank of the true polymorphism
# rule rank_true_median:
#     input:
#         filein='{macro}/{micro}/{seed}/isweep.ranks.median.tsv.gz',
#     output:
#         fileout='{macro}/{micro}/{seed}/isweep.rank.true.median.tsv',
#     params:
#         scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
#         loc=str(config['FIXED']['SIMULATE']['LOC']),
#     shell:
#         """
#         python {params.scripts}/rank-true.py \
#             {input.filein} \
#             {output.fileout} \
#             {params.loc} \
#             1
#         """

# ##### everything below this is for relate, clues #####
#
# ### first, the fast way, converting tskit to relate ###
#
# # put stuff here
#
# rule clues:
#     input:
#         anc='{macro}/{micro}/{seed}/relate_resample.anc', # ?
#         mut='{macro}/{micro}/{seed}/relate_resample.mut', # ?
#         coa='{macro}/{micro}/{seed}/relate_resample.coal',
#     output:
#         anc='{macro}/{micro}/{seed}/clues.epochs.npy',
#         mut='{macro}/{micro}/{seed}/clues.freqs.npy',
#         coa='{macro}/{micro}/{seed}/clues.post.npy',
#     params:
#         clues=str(config['CHANGE']['SOFTWARE']['CLUES']),
#         mu=str(config['FIXED']['SIMULATE']['MU']),
#         loc=str(config['FIXED']['SIMULATE']['LOC']),
#     shell:
#         """
#         python {params.clues}/inference.py \
#             --times relate_resample \
#             --coal relate_resample \
#             --out clues
#         """
#
# ### second, the slow way, running relate ###
#
# # downsample, subset to region near causal mutation
# rule relate_vcf:
#     input:
#         vcfin='{macro}/{micro}/{seed}/small.chr1.vcf.gz',
#     output:
#         vcfout='{macro}/{micro}/{seed}/relate.chr1.vcf.gz',
#     params:
#         left=str(config['FIXED']['SIMULATE']['RELATELEFT']),
#         right=str(config['FIXED']['SIMULATE']['RELATERIGHT']),
#         size=str(config['FIXED']['SIMULATE']['RELATESIZE']),
#     shell:
#         """
#         bcftools query -l {input.vcfin} | \
#             head -n {params.size} \
#             > {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/relate.keep
#         tabix -fp vcf {input.vcfin}
#         bcftools view \
#             -S {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/relate.keep
#             -r 1:{params.left}-{params.right} \
#             -Oz -o {output.vcfout} \
#             {input.vcfin}
#         """
#
# # set up input files from vcf
# rule relate_input:
#     input:
#         vcfin='{macro}/{micro}/{seed}/relate.chr1.vcf.gz',
#     output:
#         hapsout='{macro}/{micro}/{seed}/relate.haps',
#         sampout='{macro}/{micro}/{seed}/relate.sample',
#     params:
#         relate=str(config['CHANGE']['SOFTWARE']['RELATE']),
#     shell:
#         """
#         {params.relate}/bin/RelateFileFormats \
#             --mode ConvertFromVcf \
#             --haps {output.hapsout} \
#             --samples {output.sampout} \
#             --input {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/relate.chr1
#         """
#
# # run relate for first time
# rule relate_run:
#     input:
#         hapsout='{macro}/{micro}/{seed}/relate.haps',
#         sampout='{macro}/{micro}/{seed}/relate.sample',
#     output:
#         anc='{macro}/{micro}/{seed}/relate.anc',
#         mut='{macro}/{micro}/{seed}/relate.mut',
#     params:
#         relate=str(config['CHANGE']['SOFTWARE']['RELATE']),
#         mu=str(config['FIXED']['SIMULATE']['MU']),
#         size=str(config['CHANGE']['SIMULATE']['HAPSIZE']),
#     shell:
#         """
#         {params.relate}/bin/Relate \
#             --mode All \
#             --haps {output.hapsout} \
#             --samples {output.sampout} \
#             -m {params.mu} \
#             -N {params.size}\
#             -o {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/relate
#             --map {wildcards.macro}/uniform.map
#         """
#
# # create the population labels files
# rule relate_label:
#     input:
#         hapsout='{macro}/{micro}/{seed}/relate.haps',
#         sampout='{macro}/{micro}/{seed}/relate.sample',
#     output:
#         lbl='{macro}/{micro}/{seed}/relate.poplabels',
#     params:
#         scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
#     shell:
#         """
#         python {params.scripts}/poplabels.py \
#             {input.sampout} \
#             {output.lbl}
#         """
#
# # estimate effective size changes
# rule relate_pop:
#     input:
#         anc='{macro}/{micro}/{seed}/relate.anc',
#         mut='{macro}/{micro}/{seed}/relate.mut',
#         lbl='{macro}/{micro}/{seed}/relate.poplabels',
#     output:
#         anc='{macro}/{micro}/{seed}/relate_popsize.anc',
#         mut='{macro}/{micro}/{seed}/relate_popsize.mut',
#         coa='{macro}/{micro}/{seed}/relate_popsize.coal',
#     params:
#         relate=str(config['CHANGE']['SOFTWARE']['RELATE']),
#         mu=str(config['FIXED']['SIMULATE']['MU']),
#     shell:
#         """
#         {params.relate}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
#             -i {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/relate
#             -m {params.mu} \
#             --poplabels {input.lbl} \
#             -o {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/relate_popsize
#         """
#
# # resample the branch lengths
# rule relate_branch:
#     input:
#         anc='{macro}/{micro}/{seed}/relate_popsize.anc',
#         mut='{macro}/{micro}/{seed}/relate_popsize.mut',
#         coa='{macro}/{micro}/{seed}/relate_popsize.coal',
#     output:
#         anc='{macro}/{micro}/{seed}/relate_resample.anc',
#         mut='{macro}/{micro}/{seed}/relate_resample.mut',
#         coa='{macro}/{micro}/{seed}/relate_resample.coal',
#     params:
#         relate=str(config['CHANGE']['SOFTWARE']['RELATE']),
#         mu=str(config['FIXED']['SIMULATE']['MU']),
#         loc=str(config['FIXED']['SIMULATE']['LOC']),
#     shell:
#         """
#         {params.relate}/scripts/SampleBranchLengths/SampleBranchLengths.sh \
#             -i {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/relate_popsize \
#             -o {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/relate_resample \
#             -m {params.mu} \
#             --coal {input.coa} \
#             --format -b \
#             --num_samples 100 \
#             --first_bp {params.loc} --last_bp {params.loc}
#         """
#
# rule clues:
#     input:
#         anc='{macro}/{micro}/{seed}/relate_resample.anc', # ?
#         mut='{macro}/{micro}/{seed}/relate_resample.mut', # ?
#         coa='{macro}/{micro}/{seed}/relate_resample.coal',
#     output:
#         anc='{macro}/{micro}/{seed}/clues.epochs.npy',
#         mut='{macro}/{micro}/{seed}/clues.freqs.npy',
#         coa='{macro}/{micro}/{seed}/clues.post.npy',
#     params:
#         clues=str(config['CHANGE']['SOFTWARE']['CLUES']),
#         mu=str(config['FIXED']['SIMULATE']['MU']),
#         loc=str(config['FIXED']['SIMULATE']['LOC']),
#     shell:
#         """
#         python {params.clues}/inference.py \
#             --times relate_resample \
#             --coal relate_resample \
#             --out clues
#         """

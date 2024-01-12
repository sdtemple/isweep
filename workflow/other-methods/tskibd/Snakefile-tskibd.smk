# isweep in simulation study w/ true ibd
# Seth D. Temple, sdtemple@uw.edu
# October 2, 2023 added selscan

# setup macro folder
import os
macro=str(config['CHANGE']['FOLDERS']['MACRO'])
if not os.path.exists(macro):
	os.mkdir(macro)

# load in experiments, set up micro, simname folders
import pandas as pd
micro=str(config['CHANGE']["FOLDERS"]["MICRO"])
sims = pd.read_csv(micro, sep='\t', header=0)
J = sims.shape[0]
for j in range(J):
	row = sims.loc[j,]
	if not os.path.exists(macro+'/'+str(row.MICROEXP)):
		os.mkdir(macro+'/'+str(row.MICROEXP))
	if not os.path.exists(macro+'/'+str(row.MICROEXP)+'/'+str(row.SIMNAME)):
		os.mkdir(macro+'/'+str(row.MICROEXP)+'/'+str(row.SIMNAME))
sims['FOLDER'] = [macro + '/' + sims.loc[j].MICROEXP + '/' + str(sims.loc[j].SIMNAME) for j in range(J)]
sims['FILE'] = sims['FOLDER'] + '/results.hap.tsv' # after running isweep
sims['EXISTS'] = sims['FILE'].apply(os.path.isfile)
sims=sims[sims['EXISTS']]
sims = sims.set_index("SIMNAME", drop=False)

rule all:
	input:
        	# tskibd
		[f"{sim.FOLDER}/results.tskibd.tsv".replace(" ","") for sim in sims.itertuples()],
		[f"{sim.FOLDER}/tskibd.been.zipped".replace(" ","") for sim in sims.itertuples()],

# unzip a tree that had tszip applied
rule unzip_tree:
    input:
        tsz='{macro}/{micro}/{seed}/slimulation.trees.tsz',
    output:
        trees='{macro}/{micro}/{seed}/slimulation.trees',
    shell:
        """
        tsunzip -k {input.tsz}
        """

rule recap_tree:
    input:
        trees='{macro}/{micro}/{seed}/slimulation.trees',
    output:
        trees='{macro}/{micro}/{seed}/recap.trees',
    params:
        terminalscripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        rho=str(config['FIXED']['SIMULATE']['RHO']),
        ancNe=str(config['CHANGE']['SIMULATE']['ancNe']),
    shell:
        """
        python {params.terminalscripts}/recap.py \
            {input.trees} \
            {output.trees} \
            {params.rho} \
            {params.ancNe}
        """

# # not implemented
# # subset the trees file
# rule cut_tree:
#     input:
#         trees='{macro}/{micro}/{seed}/slimulation.trees',
#     output:
#         trees='{macro}/{micro}/{seed}/less.trees',
#     params:
#         terminalscripts='',
#         left='',
#         right='',
#     shell:
#         """
#         python {params.terminalscripts}/subset-interval-trees.py \
#             {input.trees} \ 
#             {output.trees} \
#             {params.left} \
#             {params.right}
#         """


# execute the main program
rule tskibd_run:
    input:
        trees='{macro}/{micro}/{seed}/recap.trees',
    output:
        ibd='{macro}/{micro}/{seed}/1.ibd',
    params:
        thresh=str(config['FIXED']['SIMULATE']['CUT']),
        binfolder=str(config['CHANGE']['FOLDERS']['TSKIBDFOLDER']),
        datfolder='{macro}/{micro}/{seed}',
    shell:
        """
        cd {params.datfolder}
        {params.binfolder}/build/tskibd \
            1 \
            1000000 \
            10000 \
            {params.thresh} \
            recap.trees      
        """

rule tskibd_gzip:
    input:
        ibd='{macro}/{micro}/{seed}/1.ibd',
    output:
        ibd='{macro}/{micro}/{seed}/1.ibd.gz',
    params:
        thresh=str(config['FIXED']['SIMULATE']['CUT']),
        binfolder=str(config['CHANGE']['FOLDERS']['TSKIBDFOLDER']),
        datfolder='{macro}/{micro}/{seed}',
    shell:
        """
        gzip -c {input.ibd} > {output.ibd}       
        """

# filter some lines
rule tskibd_filter:
    input:
        ibd='{macro}/{micro}/{seed}/1.ibd',
        ibdgz='{macro}/{micro}/{seed}/1.ibd.gz',
    output:
        ibd='{macro}/{micro}/{seed}/third.tskibd.ibd.gz',
    params:
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        loc=str(config['FIXED']['SIMULATE']['LOC']),
    shell:
        """
        cat {input.ibd} | \
            java -jar {params.soft}/filter-lines.jar \
            "I" 3 0.00 {params.loc} | \
            java -jar {params.soft}/filter-lines.jar \
            "I" 4 {params.loc} 10000000000 | \
            gzip > {output.ibd}
        rm {input.ibd}
        """

# requires isweep
rule tskibd_infer:
    input:
        long='{macro}/{micro}/{seed}/third.tskibd.ibd.gz',
        freq='{macro}/{micro}/{seed}/third.best.hap.txt',
    output:
        fileout='{macro}/{micro}/{seed}/results.tskibd.tsv',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        nboot=str(config['FIXED']['ISWEEP']['NBOOT']),
        cm=str(config['FIXED']['ISWEEP']['MOMCUTOFF']),
        n=str(config['FIXED']['SIMULATE']['SAMPSIZE']),
        ploidy=str(config['FIXED']['HAPIBD']['PLOIDY']),
        effdemo=str(config['CHANGE']['SIMULATE']['iNe']),
    shell:
        """
        ibdest=$(zcat {input.long} | wc -l)
        freqest=$(python {params.scripts}/lines.py {input.freq} 2 2)
        python {params.scripts}/estimate.py \
            {output.fileout} \
            ${{ibdest}} \
            ${{freqest}} \
            {params.nboot} \
            {params.cm} \
            {params.n} \
            {params.effdemo} \
            {params.ploidy}
        """

rule zip_tree:
    input:
        finished='{macro}/{micro}/{seed}/results.tskibd.tsv',
    output:
        tsz='{macro}/{micro}/{seed}/slimulation.trees.tsz',
        rtsz='{macro}/{micro}/{seed}/recap.trees.tsz',
        tszipped='{macro}/{micro}/{seed}/tskibd.been.zipped'
    params:
        trees='{macro}/{micro}/{seed}/slimulation.trees',
        trees2='{macro}/{micro}/{seed}/recap.trees',
    shell:
        """
        tszip {params.trees}
        tszip {params.trees2}
        touch {output.tszipped}  
        """

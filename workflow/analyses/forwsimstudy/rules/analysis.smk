wildcard_constraints:
	SIMNAME = '\w+',

rule community:
    input:
        ibdin='{macro}/{micro}/small.{method}.ibd.gz',
        vcfin='{macro}/{micro}/small.{method}.vcf.gz',
    params:
        header='0',
    output:
        fileout='{macro}/{micro}/short.{method}.full.ranks.tsv.gz',
    script:
        '../../../snakescripts/isweep-rank.py'

# rule kmeans:
# 	input:
# 		filein='{macro}/{micro}/short.{method}.full.ranks.tsv.gz',
# 	output:
# 		fileout='{macro}/{micro}/short.{method}.kmeans.ranks.tsv',
# 		fignout='{macro}/{micro}/short.{method}.ranks.png',
# 	script:
# 		'../../../snakescripts/isweep-cluster.py'

rule analysis:
	input:
		ibdlong='{macro}/{micro}/long.{method}.ibd.gz',
		ibdcomm='{macro}/{micro}/short.{method}.full.ranks.tsv.gz',
		freq='{macro}/{micro}/slimulation.freq',
	output:
		fileout='{macro}/{micro}/{method}.results.tsv',
	params:
		header='0',
	script:
		'../../../snakescripts/isweep-analysis.py'

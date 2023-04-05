wildcard_constraints:
	SIMNAME = '\w+',

# these prefixes should match that in rule all in Snakefile
# pipeline will not run otherwise
rule remove:
	input:
		'{macro}/{micro}/hapibd.results.tsv',
		# '{macro}/{micro}/ibdends.true.results.tsv',
	output:
		touched='{macro}/{micro}/removed',
	shell:
		"""
		rm -f {macro}/{micro}/large.vcf.gz
		rm -f {macro}/{micro}/slimulation.bcf.gz.csi
		rm -f {macro}/{micro}/slimulation.bcf.gz
		rm -f {macro}/{micro}/slimulation.vcf.gz
		rm -f {macro}/{micro}/large.vcf.bgz.tbi
		rm -f {macro}/{micro}/large.vcf.bgz
		rm -f {macro}/{micro}/short.hapibd.ibd.gz
		rm -f {macro}/{micro}/short.ibdends.ibd.gz
		rm -f {macro}/{micro}/long.hapibd.ibd.gz
		rm -f {macro}/{micro}/long.ibdends.ibd.gz
		rm -f {macro}/{micro}/large.ibdends.p1.ibd.gz
		rm -f {macro}/{micro}/large.ibdends.p1.hbd.gz
		rm -f {macro}/{micro}/large.hapibd.cand.ibd.gz
		rm -f {macro}/{micro}/large.hapibd.cand.hbd.gz
		rm -f {macro}/{micro}/firstpass.ibdends.ibd.gz
		rm -f {macro}/{micro}/firstpass.ibdends.ibd.gz
		rm -f {macro}/{micro}/large.hapibd.ibd.gz
		rm -f {macro}/{micro}/large.ibdends.ibd.gz
		rm -f {macro}/{micro}/large.hapibd.hbd.gz
		rm -f {macro}/{micro}/large.hapibd.vcf.bgz.tbi
		rm -f {macro}/{micro}/large.hapibd.vcf.bgz
		rm -f {macro}/{micro}/large.ibdends.vcf.bgz.tbi
		rm -f {macro}/{micro}/large.ibdends.vcf.bgz
		touch {output.touched}
		"""

rule record:
	input:
		'{macro}/{micro}/removed',
	output:
		yaml='{macro}/{micro}/arguments.yaml'
	shell:
		'cp ../config/simstudies.yaml {output.yaml}'

rule tszip:
	input:
		trees='{macro}/{micro}/slimulation.trees',
		yaml='{macro}/{micro}/arguments.yaml',
	output:
		tsz='{macro}/{micro}/slimulation.trees.tsz'
	shell:
		'tszip {input.trees}'

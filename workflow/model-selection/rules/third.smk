### estimating selection coefficients
### estimating allele frequencies

# some inputs, string managements, count sample size
cohort=str(config['CHANGE']['FOLDERS']['STUDY'])
samplesize=0
with open(cohort+'/subsample.txt','r') as f:
    for line in f:
        samplesize+=1
ploidy=2
# ploidy=int(float(str(config['FIXED']['HAPIBD']['PLOIDY'])))
samplesize_ploidy=samplesize*ploidy

import os
import pandas as pd
macro=str(config['CHANGE']['FOLDERS']['STUDY'])
micro=str(config['CHANGE']["ISWEEP"]["ROI"])
sims = pd.read_csv(macro+'/'+micro, sep='\t', header=0)
J = sims.shape[0]
for j in range(J):
	row = sims.iloc[j]
	if not os.path.exists(macro+'/'+str(row.NAME)):
		os.mkdir(macro+'/'+str(row.NAME))
sims['FOLDER'] = [(macro +'/'+str(sims.iloc[j].NAME)).strip() for j in range(J)]

# extend Ne(t)
rule third_Ne:
    input:
        shortNe='{cohort}/'+str(config['CHANGE']['ISWEEP']['NE']),
    output:
        longNe='{cohort}/extended.ne',
    params:
        lastne=str(config['FIXED']['ISWEEP']['LASTNE']),
    shell:
        """
        python -c "from isweep import extend_Ne; extend_Ne('{input.shortNe}',float('{params.lastne}'),'{output.longNe}')"
        """

##### haplotype analysis #####

rule third_hap:
    input:
        rankin='{cohort}/{hit}/second.ranks.tsv.gz',
        outlied='{cohort}/{hit}/second.outliers.txt',
    output:
        lociout='{cohort}/{hit}/third.best.hap.txt',
        happng='{cohort}/{hit}/third.hap.png',
        snppng='{cohort}/{hit}/third.snp.png',
    params:
        windowsize=str(config['FIXED']['ISWEEP']['HAPSIZE']),
        windowstep=str(config['FIXED']['ISWEEP']['HAPSTEP']),
        freqsize=str(config['FIXED']['ISWEEP']['FREQSIZE']),
        freqstep=str(config['FIXED']['ISWEEP']['FREQSTEP']),
        numsnp=str(config['FIXED']['ISWEEP']['NUMSNP']),
        lowbnd=str(config['FIXED']['ISWEEP']['MINAAF']),
    shell:
        """
        python ../../scripts/model/haplotypes.py \
            --input_snp_file {input.rankin} \
            --output_folder {wildcards.cohort}/{wildcards.hit} \
            --window_index 0 \
            --freq_index 1 \
            --score_index -1 \
            --freq_size {params.freqsize} \
            --freq_step {params.freqstep} \
            --window_size {params.windowsize} \
            --window_step {params.windowstep} \
            --num_snp {params.numsnp} \
            --lowest_freq {params.lowbnd}
        """


rule third_hap_ibd:
    input:
        best='{cohort}/{hit}/third.best.hap.txt',
        locus='{cohort}/{hit}/locus.txt',
    output:
        ibd='{cohort}/{hit}/third.hap.ibd.gz',
    params:
        ibdfolder='{cohort}/ibdsegs/ibdends/mle',
    shell:
        """
        chr=$(python ../../scripts/utilities/lines.py {input.locus} 2 2)
        thecenter=$(python ../../scripts/utilities/lines.py {input.best} 1 2)
        ibd={params.ibdfolder}/chr${{chr}}.ibd.gz
        python ../../scripts/utilities/filter-lines.py \
            --input_file $ibd \
            --output_file {wildcards.cohort}/{wildcards.hit}/intermediate.ibd.gz \
            --column_index 6 \
            --upper_bound $thecenter \
            --complement 0
        python ../../scripts/utilities/filter-lines.py \
            --input_file {wildcards.cohort}/{wildcards.hit}/intermediate.ibd.gz \
            --output_file {output.ibd} \
            --column_index 7 \
            --lower_bound $thecenter \
            --upper_bound 10000000000 \
            --complement 0
        rm {wildcards.cohort}/{wildcards.hit}/intermediate.ibd.gz
        """


rule third_hap_infer_norm:
    input:
        long='{cohort}/{hit}/third.hap.ibd.gz',
        freq='{cohort}/{hit}/third.best.hap.txt',
        loci='{cohort}/{hit}/locus.txt',
        longNe='{cohort}/extended.ne',
    output:
        fileout='{cohort}/{hit}/results.hap.norm.tsv',
    params:
        nboot=str(config['FIXED']['ISWEEP']['NBOOTNORM']),
        cm=str(config['FIXED']['ISWEEP']['MLECUTOFF']),
        n=str(samplesize),
        ploidy=str(2),
        # ploidy=str(config['FIXED']['HAPIBD']['PLOIDY'])
    shell:
        """
        ibdest=$(zcat {input.long} | wc -l)
        freqest=$(python ../../scripts/utilities/lines.py {input.freq} 2 2)
        model=$(python ../../scripts/utilities/lines.py {input.loci} 6 2)
        alpha=$(python ../../scripts/utilities/lines.py {input.loci} 7 2)
        python ../../scripts/model/estimate-norm.py \
            --output_file {output.fileout} \
            --ibd_count ${{ibdest}} \
            --p_est ${{freqest}} \
            --num_bootstraps {params.nboot} \
            --ibd_cutoff {params.cm} \
            --sample_size {params.n} \
            --Ne_est {input.longNe} \
            --model ${{model}} \
            --alpha ${{alpha}} \
            --ploidy {params.ploidy} \
        """

rule third_hap_infer_perc:
    input:
        long='{cohort}/{hit}/third.hap.ibd.gz',
        freq='{cohort}/{hit}/third.best.hap.txt',
        loci='{cohort}/{hit}/locus.txt',
        longNe='{cohort}/extended.ne',
    output:
        fileout='{cohort}/{hit}/results.hap.perc.tsv',
    params:
        nboot=str(config['FIXED']['ISWEEP']['NBOOTPERC']),
        cm=str(config['FIXED']['ISWEEP']['MLECUTOFF']),
        n=str(samplesize),
        ploidy=str(2),
        # ploidy=str(config['FIXED']['HAPIBD']['PLOIDY'])
    shell:
        """
        ibdest=$(zcat {input.long} | wc -l)
        freqest=$(python ../../scripts/utilities/lines.py {input.freq} 2 2)
        model=$(python ../../scripts/utilities/lines.py {input.loci} 6 2)
        alpha=$(python ../../scripts/utilities/lines.py {input.loci} 7 2)
        python ../../scripts/model/estimate-perc.py \
            --output_file {output.fileout} \
            --ibd_count ${{ibdest}} \
            --p_est ${{freqest}} \
            --num_bootstraps {params.nboot} \
            --ibd_cutoff {params.cm} \
            --sample_size {params.n} \
            --Ne_est {input.longNe} \
            --model ${{model}} \
            --alpha ${{alpha}} \
            --ploidy {params.ploidy} \
        """

rule summary_hap_norm:
    input:
        [(macro +'/'+str(sims.iloc[j].NAME)).strip()+'/results.hap.norm.tsv' for j in range(J)],
        [(macro +'/'+str(sims.iloc[j].NAME)).strip()+'/ibd.gini.tsv' for j in range(J)],
    output:
        fileout=macro+'/summary.hap.norm.tsv',
    params:
        study=macro,
        roi=str(config['CHANGE']['ISWEEP']['ROI']),
        typ='hap',
        unc='norm',
    shell:
        """
        python ../../scripts/model/summary-table.py \
            --output_file {output.fileout} \
            --input_folder {params.study} \
            --input_roi_file {params.roi} \
            --file_type {params.typ} \
            --uncertainty_type {params.unc} \
            --input_gini_file ibd.gini.tsv \
        """

rule summary_hap_perc:
    input:
        [(macro +'/'+str(sims.iloc[j].NAME)).strip()+'/results.hap.perc.tsv' for j in range(J)],
        [(macro +'/'+str(sims.iloc[j].NAME)).strip()+'/ibd.gini.tsv' for j in range(J)],
    output:
        fileout=macro+'/summary.hap.perc.tsv',
    params:
        study=macro,
        roi=str(config['CHANGE']['ISWEEP']['ROI']),
        typ='hap',
        unc='perc',
    shell:
        """
        python ../../scripts/model/summary-table.py \
            --output_file {output.fileout} \
            --input_folder {params.study} \
            --input_roi_file {params.roi} \
            --file_type {params.typ} \
            --uncertainty_type {params.unc} \
            --input_gini_file ibd.gini.tsv \
        """

##### snp analysis #####

rule third_snp:
    input:
        rankin='{cohort}/{hit}/second.ranks.tsv.gz',
        outlied='{cohort}/{hit}/second.outliers.txt',
    output:
        lociout='{cohort}/{hit}/third.best.snp.txt',
    params:
        lowbnd=str(config['FIXED']['ISWEEP']['MINAAF']),
    shell:
        """
        python ../../scripts/model/snp.py \
            --input_snp_file {input.rankin} \
            --output_file {output.lociout} \
            --lowest_freq {params.lowbnd}
        """

rule third_snp_ibd:
    input:
        best='{cohort}/{hit}/third.best.snp.txt',
        locus='{cohort}/{hit}/locus.txt',
    output:
        ibd='{cohort}/{hit}/third.snp.ibd.gz',
    params:
        ibdfolder='{cohort}/ibdsegs/ibdends/mle',
    shell:
        """
        chr=$(python ../../scripts/utilities/lines.py {input.locus} 2 2)
        thecenter=$(python ../../scripts/utilities/lines.py {input.best} 1 2)
        ibd={params.ibdfolder}/chr${{chr}}.ibd.gz
        python ../../scripts/utilities/filter-lines.py \
            --input_file $ibd \
            --output_file {wildcards.cohort}/{wildcards.hit}/intermediate2.ibd.gz \
            --column_index 6 \
            --upper_bound $thecenter \
            --complement 0
        python ../../scripts/utilities/filter-lines.py \
            --input_file {wildcards.cohort}/{wildcards.hit}/intermediate2.ibd.gz \
            --output_file {output.ibd} \
            --column_index 7 \
            --lower_bound $thecenter \
            --upper_bound 10000000000 \
            --complement 0
        rm {wildcards.cohort}/{wildcards.hit}/intermediate2.ibd.gz
        """


rule third_snp_infer_norm:
    input:
        long='{cohort}/{hit}/third.snp.ibd.gz',
        freq='{cohort}/{hit}/third.best.snp.txt',
        loci='{cohort}/{hit}/locus.txt',
        longNe='{cohort}/extended.ne',
    output:
        fileout='{cohort}/{hit}/results.snp.norm.tsv',
    params:
        nboot=str(config['FIXED']['ISWEEP']['NBOOTNORM']),
        cm=str(config['FIXED']['ISWEEP']['MLECUTOFF']),
        n=str(samplesize),
        ploidy=str(2),
        # ploidy=str(config['FIXED']['HAPIBD']['PLOIDY'])
    shell:
        """
        ibdest=$(zcat {input.long} | wc -l)
        freqest=$(python ../../scripts/utilities/lines.py {input.freq} 2 2)
        model=$(python ../../scripts/utilities/lines.py {input.loci} 6 2)
        alpha=$(python ../../scripts/utilities/lines.py {input.loci} 7 2)
        python ../../scripts/model/estimate-norm.py \
            --output_file {output.fileout} \
            --ibd_count ${{ibdest}} \
            --p_est ${{freqest}} \
            --num_bootstraps {params.nboot} \
            --ibd_cutoff {params.cm} \
            --sample_size {params.n} \
            --Ne_est {input.longNe} \
            --model ${{model}} \
            --alpha ${{alpha}} \
            --ploidy {params.ploidy} \
        """

rule third_snp_infer_perc:
    input:
        long='{cohort}/{hit}/third.snp.ibd.gz',
        freq='{cohort}/{hit}/third.best.snp.txt',
        loci='{cohort}/{hit}/locus.txt',
        longNe='{cohort}/extended.ne',
    output:
        fileout='{cohort}/{hit}/results.snp.perc.tsv',
    params:
        nboot=str(config['FIXED']['ISWEEP']['NBOOTPERC']),
        cm=str(config['FIXED']['ISWEEP']['MLECUTOFF']),
        n=str(samplesize),
        ploidy=str(2),
        # ploidy=str(config['FIXED']['HAPIBD']['PLOIDY'])
    shell:
        """
        ibdest=$(zcat {input.long} | wc -l)
        freqest=$(python ../../scripts/utilities/lines.py {input.freq} 2 2)
        model=$(python ../../scripts/utilities/lines.py {input.loci} 6 2)
        alpha=$(python ../../scripts/utilities/lines.py {input.loci} 7 2)
        python ../../scripts/model/estimate-perc.py \
            --output_file {output.fileout} \
            --ibd_count ${{ibdest}} \
            --p_est ${{freqest}} \
            --num_bootstraps {params.nboot} \
            --ibd_cutoff {params.cm} \
            --sample_size {params.n} \
            --Ne_est {input.longNe} \
            --model ${{model}} \
            --alpha ${{alpha}} \
            --ploidy {params.ploidy} \
        """

rule summary_snp_norm:
    input:
        [(macro +'/'+str(sims.iloc[j].NAME)).strip()+'/results.snp.norm.tsv' for j in range(J)],
        [(macro +'/'+str(sims.iloc[j].NAME)).strip()+'/ibd.gini.tsv' for j in range(J)],
    output:
        fileout=macro+'/summary.snp.norm.tsv',
    params:
        study=macro,
        roi=str(config['CHANGE']['ISWEEP']['ROI']),
        typ='snp',
        unc='norm',
    shell:
        """
        python ../../scripts/model/summary-table.py \
            --output_file {output.fileout} \
            --input_folder {params.study} \
            --input_roi_file {params.roi} \
            --file_type {params.typ} \
            --uncertainty_type {params.unc} \
            --input_gini_file ibd.gini.tsv \
        """

rule summary_snp_perc:
    input:
        [(macro +'/'+str(sims.iloc[j].NAME)).strip()+'/results.snp.perc.tsv' for j in range(J)],
        [(macro +'/'+str(sims.iloc[j].NAME)).strip()+'/ibd.gini.tsv' for j in range(J)],
    output:
        fileout=macro+'/summary.snp.perc.tsv',
    params:
        study=macro,
        roi=str(config['CHANGE']['ISWEEP']['ROI']),
        typ='snp',
        unc='perc',
    shell:
        """
        python ../../scripts/model/summary-table.py \
            --output_file {output.fileout} \
            --input_folder {params.study} \
            --input_roi_file {params.roi} \
            --file_type {params.typ} \
            --uncertainty_type {params.unc} \
            --input_gini_file ibd.gini.tsv \
        """

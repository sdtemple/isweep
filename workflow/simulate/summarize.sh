# summarize snakemake simulate workflow
# seth temple
# july 18, 2023

macro=$1
rm -r $macro/concat/
rm -r $macro/forplots/
rm -r $macro/fortables/
mkdir -p $macro/concat/
mkdir -p $macro/forplots/
mkdir -p $macro/fortables/

micro=MIDsEARLYt
left=1
right=100
# write files
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/first.pos.txt | cut -f 2 >> $macro/concat/$micro.second.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/first.mode.txt | cut -f 2 >> $macro/concat/$micro.first.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do head -n 1 $macro/$micro/$j/third.best.hap.txt | cut -f 2 >> $macro/concat/$micro.hap.third.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do head -n 1 $macro/$micro/$j/third.best.snp.txt | cut -f 2 >> $macro/concat/$micro.snp.third.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 2-4 >> $macro/concat/$micro.hap.sel.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 2-4 >> $macro/concat/$micro.snp.sel.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 1 >> $macro/concat/$micro.hap.freq.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 1 >> $macro/concat/$micro.snp.freq.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 6-8 >> $macro/concat/$micro.hap.time5.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 6-8 >> $macro/concat/$micro.snp.time5.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 9-11 >> $macro/concat/$micro.hap.time2.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 9-11 >> $macro/concat/$micro.snp.time2.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 12-14 >> $macro/concat/$micro.hap.denovo.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 12-14 >> $macro/concat/$micro.snp.denovo.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/mathieson.txt >> $macro/concat/$micro.mathieson.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/slimulation.freq >> $macro/concat/$micro.slim.freq.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/clues.v1.selcoef >> $macro/concat/$micro.clues.txt ; done ;
# paste files
paste $macro/concat/$micro.second.mode.txt $macro/concat/$micro.first.mode.txt $macro/concat/$micro.hap.third.mode.txt $macro/concat/$micro.snp.third.mode.txt > $macro/forplots/$micro.modes.tsv
paste $macro/concat/$micro.hap.freq.txt $macro/concat/$micro.snp.freq.txt $macro/concat/$micro.slim.freq.txt > $macro/forplots/$micro.freq.tsv
paste $macro/concat/$micro.hap.sel.tsv $macro/concat/$micro.snp.sel.tsv $macro/concat/$micro.mathieson.txt $macro/concat/$micro.clues.txt | cut -f 1,4,7,10 > $macro/forplots/$micro.sel.plots.tsv
paste $macro/concat/$micro.hap.sel.tsv $macro/concat/$micro.snp.sel.tsv $macro/concat/$micro.mathieson.txt > $macro/fortables/$micro.sel.tables.tsv
paste $macro/concat/$micro.hap.time5.tsv $macro/concat/$micro.hap.time2.tsv $macro/concat/$micro.hap.denovo.tsv > $macro/fortables/$micro.hap.time.tables.tsv
paste $macro/concat/$micro.snp.time5.tsv $macro/concat/$micro.snp.time2.tsv $macro/concat/$micro.snp.denovo.tsv > $macro/fortables/$micro.snp.time.tables.tsv


micro=MIDsLATEt
left=1
right=100
# write files
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/first.pos.txt | cut -f 2 >> $macro/concat/$micro.second.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/first.mode.txt | cut -f 2 >> $macro/concat/$micro.first.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do head -n 1 $macro/$micro/$j/third.best.hap.txt | cut -f 2 >> $macro/concat/$micro.hap.third.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do head -n 1 $macro/$micro/$j/third.best.snp.txt | cut -f 2 >> $macro/concat/$micro.snp.third.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 2-4 >> $macro/concat/$micro.hap.sel.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 2-4 >> $macro/concat/$micro.snp.sel.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 1 >> $macro/concat/$micro.hap.freq.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 1 >> $macro/concat/$micro.snp.freq.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 6-8 >> $macro/concat/$micro.hap.time5.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 6-8 >> $macro/concat/$micro.snp.time5.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 9-11 >> $macro/concat/$micro.hap.time2.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 9-11 >> $macro/concat/$micro.snp.time2.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 12-14 >> $macro/concat/$micro.hap.denovo.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 12-14 >> $macro/concat/$micro.snp.denovo.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/mathieson.txt >> $macro/concat/$micro.mathieson.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/slimulation.freq >> $macro/concat/$micro.slim.freq.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/clues.v1.selcoef >> $macro/concat/$micro.clues.txt ; done ;
# paste files
paste $macro/concat/$micro.second.mode.txt $macro/concat/$micro.first.mode.txt $macro/concat/$micro.hap.third.mode.txt $macro/concat/$micro.snp.third.mode.txt > $macro/forplots/$micro.modes.tsv
paste $macro/concat/$micro.hap.freq.txt $macro/concat/$micro.snp.freq.txt $macro/concat/$micro.slim.freq.txt > $macro/forplots/$micro.freq.tsv
paste $macro/concat/$micro.hap.sel.tsv $macro/concat/$micro.snp.sel.tsv $macro/concat/$micro.mathieson.txt $macro/concat/$micro.clues.txt | cut -f 1,4,7,10 > $macro/forplots/$micro.sel.plots.tsv
paste $macro/concat/$micro.hap.sel.tsv $macro/concat/$micro.snp.sel.tsv $macro/concat/$micro.mathieson.txt > $macro/fortables/$micro.sel.tables.tsv
paste $macro/concat/$micro.hap.time5.tsv $macro/concat/$micro.hap.time2.tsv $macro/concat/$micro.hap.denovo.tsv > $macro/fortables/$micro.hap.time.tables.tsv
paste $macro/concat/$micro.snp.time5.tsv $macro/concat/$micro.snp.time2.tsv $macro/concat/$micro.snp.denovo.tsv > $macro/fortables/$micro.snp.time.tables.tsv

micro=TINYsEARLYt
left=1
right=100
# write files
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/first.pos.txt | cut -f 2 >> $macro/concat/$micro.second.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/first.mode.txt | cut -f 2 >> $macro/concat/$micro.first.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do head -n 1 $macro/$micro/$j/third.best.hap.txt | cut -f 2 >> $macro/concat/$micro.hap.third.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do head -n 1 $macro/$micro/$j/third.best.snp.txt | cut -f 2 >> $macro/concat/$micro.snp.third.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 2-4 >> $macro/concat/$micro.hap.sel.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 2-4 >> $macro/concat/$micro.snp.sel.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 1 >> $macro/concat/$micro.hap.freq.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 1 >> $macro/concat/$micro.snp.freq.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 6-8 >> $macro/concat/$micro.hap.time5.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 6-8 >> $macro/concat/$micro.snp.time5.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 9-11 >> $macro/concat/$micro.hap.time2.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 9-11 >> $macro/concat/$micro.snp.time2.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 12-14 >> $macro/concat/$micro.hap.denovo.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 12-14 >> $macro/concat/$micro.snp.denovo.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/mathieson.txt >> $macro/concat/$micro.mathieson.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/slimulation.freq >> $macro/concat/$micro.slim.freq.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/clues.v1.selcoef >> $macro/concat/$micro.clues.txt ; done ;
# paste files
paste $macro/concat/$micro.second.mode.txt $macro/concat/$micro.first.mode.txt $macro/concat/$micro.hap.third.mode.txt $macro/concat/$micro.snp.third.mode.txt > $macro/forplots/$micro.modes.tsv
paste $macro/concat/$micro.hap.freq.txt $macro/concat/$micro.snp.freq.txt $macro/concat/$micro.slim.freq.txt > $macro/forplots/$micro.freq.tsv
paste $macro/concat/$micro.hap.sel.tsv $macro/concat/$micro.snp.sel.tsv $macro/concat/$micro.mathieson.txt $macro/concat/$micro.clues.txt | cut -f 1,4,7,10 > $macro/forplots/$micro.sel.plots.tsv
paste $macro/concat/$micro.hap.sel.tsv $macro/concat/$micro.snp.sel.tsv $macro/concat/$micro.mathieson.txt > $macro/fortables/$micro.sel.tables.tsv
paste $macro/concat/$micro.hap.time5.tsv $macro/concat/$micro.hap.time2.tsv $macro/concat/$micro.hap.denovo.tsv > $macro/fortables/$micro.hap.time.tables.tsv
paste $macro/concat/$micro.snp.time5.tsv $macro/concat/$micro.snp.time2.tsv $macro/concat/$micro.snp.denovo.tsv > $macro/fortables/$micro.snp.time.tables.tsv

micro=TINYsLATEt
left=101
right=200
# write files
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/first.pos.txt | cut -f 2 >> $macro/concat/$micro.second.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/first.mode.txt | cut -f 2 >> $macro/concat/$micro.first.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do head -n 1 $macro/$micro/$j/third.best.hap.txt | cut -f 2 >> $macro/concat/$micro.hap.third.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do head -n 1 $macro/$micro/$j/third.best.snp.txt | cut -f 2 >> $macro/concat/$micro.snp.third.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 2-4 >> $macro/concat/$micro.hap.sel.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 2-4 >> $macro/concat/$micro.snp.sel.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 1 >> $macro/concat/$micro.hap.freq.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 1 >> $macro/concat/$micro.snp.freq.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 6-8 >> $macro/concat/$micro.hap.time5.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 6-8 >> $macro/concat/$micro.snp.time5.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 9-11 >> $macro/concat/$micro.hap.time2.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 9-11 >> $macro/concat/$micro.snp.time2.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 12-14 >> $macro/concat/$micro.hap.denovo.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 12-14 >> $macro/concat/$micro.snp.denovo.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/mathieson.txt >> $macro/concat/$micro.mathieson.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/slimulation.freq >> $macro/concat/$micro.slim.freq.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/clues.v1.selcoef >> $macro/concat/$micro.clues.txt ; done ;
# paste files
paste $macro/concat/$micro.second.mode.txt $macro/concat/$micro.first.mode.txt $macro/concat/$micro.hap.third.mode.txt $macro/concat/$micro.snp.third.mode.txt > $macro/forplots/$micro.modes.tsv
paste $macro/concat/$micro.hap.freq.txt $macro/concat/$micro.snp.freq.txt $macro/concat/$micro.slim.freq.txt > $macro/forplots/$micro.freq.tsv
paste $macro/concat/$micro.hap.sel.tsv $macro/concat/$micro.snp.sel.tsv $macro/concat/$micro.mathieson.txt $macro/concat/$micro.clues.txt | cut -f 1,4,7,10 > $macro/forplots/$micro.sel.plots.tsv
paste $macro/concat/$micro.hap.sel.tsv $macro/concat/$micro.snp.sel.tsv $macro/concat/$micro.mathieson.txt > $macro/fortables/$micro.sel.tables.tsv
paste $macro/concat/$micro.hap.time5.tsv $macro/concat/$micro.hap.time2.tsv $macro/concat/$micro.hap.denovo.tsv > $macro/fortables/$micro.hap.time.tables.tsv
paste $macro/concat/$micro.snp.time5.tsv $macro/concat/$micro.snp.time2.tsv $macro/concat/$micro.snp.denovo.tsv > $macro/fortables/$micro.snp.time.tables.tsv

micro=SMALLsEARLYt
left=1
right=100
# write files
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/first.pos.txt | cut -f 2 >> $macro/concat/$micro.second.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/first.mode.txt | cut -f 2 >> $macro/concat/$micro.first.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do head -n 1 $macro/$micro/$j/third.best.hap.txt | cut -f 2 >> $macro/concat/$micro.hap.third.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do head -n 1 $macro/$micro/$j/third.best.snp.txt | cut -f 2 >> $macro/concat/$micro.snp.third.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 2-4 >> $macro/concat/$micro.hap.sel.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 2-4 >> $macro/concat/$micro.snp.sel.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 1 >> $macro/concat/$micro.hap.freq.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 1 >> $macro/concat/$micro.snp.freq.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 6-8 >> $macro/concat/$micro.hap.time5.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 6-8 >> $macro/concat/$micro.snp.time5.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 9-11 >> $macro/concat/$micro.hap.time2.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 9-11 >> $macro/concat/$micro.snp.time2.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 12-14 >> $macro/concat/$micro.hap.denovo.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 12-14 >> $macro/concat/$micro.snp.denovo.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/mathieson.txt >> $macro/concat/$micro.mathieson.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/slimulation.freq >> $macro/concat/$micro.slim.freq.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/clues.v1.selcoef >> $macro/concat/$micro.clues.txt ; done ;
# paste files
paste $macro/concat/$micro.second.mode.txt $macro/concat/$micro.first.mode.txt $macro/concat/$micro.hap.third.mode.txt $macro/concat/$micro.snp.third.mode.txt > $macro/forplots/$micro.modes.tsv
paste $macro/concat/$micro.hap.freq.txt $macro/concat/$micro.snp.freq.txt $macro/concat/$micro.slim.freq.txt > $macro/forplots/$micro.freq.tsv
paste $macro/concat/$micro.hap.sel.tsv $macro/concat/$micro.snp.sel.tsv $macro/concat/$micro.mathieson.txt $macro/concat/$micro.clues.txt | cut -f 1,4,7,10 > $macro/forplots/$micro.sel.plots.tsv
paste $macro/concat/$micro.hap.sel.tsv $macro/concat/$micro.snp.sel.tsv $macro/concat/$micro.mathieson.txt > $macro/fortables/$micro.sel.tables.tsv
paste $macro/concat/$micro.hap.time5.tsv $macro/concat/$micro.hap.time2.tsv $macro/concat/$micro.hap.denovo.tsv > $macro/fortables/$micro.hap.time.tables.tsv
paste $macro/concat/$micro.snp.time5.tsv $macro/concat/$micro.snp.time2.tsv $macro/concat/$micro.snp.denovo.tsv > $macro/fortables/$micro.snp.time.tables.tsv

micro=SMALLsLATEt
left=1
right=100
# write files
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/first.pos.txt | cut -f 2 >> $macro/concat/$micro.second.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/first.mode.txt | cut -f 2 >> $macro/concat/$micro.first.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do head -n 1 $macro/$micro/$j/third.best.hap.txt | cut -f 2 >> $macro/concat/$micro.hap.third.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do head -n 1 $macro/$micro/$j/third.best.snp.txt | cut -f 2 >> $macro/concat/$micro.snp.third.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 2-4 >> $macro/concat/$micro.hap.sel.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 2-4 >> $macro/concat/$micro.snp.sel.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 1 >> $macro/concat/$micro.hap.freq.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 1 >> $macro/concat/$micro.snp.freq.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 6-8 >> $macro/concat/$micro.hap.time5.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 6-8 >> $macro/concat/$micro.snp.time5.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 9-11 >> $macro/concat/$micro.hap.time2.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 9-11 >> $macro/concat/$micro.snp.time2.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 12-14 >> $macro/concat/$micro.hap.denovo.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 12-14 >> $macro/concat/$micro.snp.denovo.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/mathieson.txt >> $macro/concat/$micro.mathieson.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/slimulation.freq >> $macro/concat/$micro.slim.freq.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/clues.v1.selcoef >> $macro/concat/$micro.clues.txt ; done ;
# paste files
paste $macro/concat/$micro.second.mode.txt $macro/concat/$micro.first.mode.txt $macro/concat/$micro.hap.third.mode.txt $macro/concat/$micro.snp.third.mode.txt > $macro/forplots/$micro.modes.tsv
paste $macro/concat/$micro.hap.freq.txt $macro/concat/$micro.snp.freq.txt $macro/concat/$micro.slim.freq.txt > $macro/forplots/$micro.freq.tsv
paste $macro/concat/$micro.hap.sel.tsv $macro/concat/$micro.snp.sel.tsv $macro/concat/$micro.mathieson.txt $macro/concat/$micro.clues.txt | cut -f 1,4,7,10 > $macro/forplots/$micro.sel.plots.tsv
paste $macro/concat/$micro.hap.sel.tsv $macro/concat/$micro.snp.sel.tsv $macro/concat/$micro.mathieson.txt > $macro/fortables/$micro.sel.tables.tsv
paste $macro/concat/$micro.hap.time5.tsv $macro/concat/$micro.hap.time2.tsv $macro/concat/$micro.hap.denovo.tsv > $macro/fortables/$micro.hap.time.tables.tsv
paste $macro/concat/$micro.snp.time5.tsv $macro/concat/$micro.snp.time2.tsv $macro/concat/$micro.snp.denovo.tsv > $macro/fortables/$micro.snp.time.tables.tsv

micro=LARGEsEARLYt
left=1
right=100
# write files
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/first.pos.txt | cut -f 2 >> $macro/concat/$micro.second.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/first.mode.txt | cut -f 2 >> $macro/concat/$micro.first.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do head -n 1 $macro/$micro/$j/third.best.hap.txt | cut -f 2 >> $macro/concat/$micro.hap.third.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do head -n 1 $macro/$micro/$j/third.best.snp.txt | cut -f 2 >> $macro/concat/$micro.snp.third.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 2-4 >> $macro/concat/$micro.hap.sel.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 2-4 >> $macro/concat/$micro.snp.sel.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 1 >> $macro/concat/$micro.hap.freq.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 1 >> $macro/concat/$micro.snp.freq.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 6-8 >> $macro/concat/$micro.hap.time5.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 6-8 >> $macro/concat/$micro.snp.time5.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 9-11 >> $macro/concat/$micro.hap.time2.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 9-11 >> $macro/concat/$micro.snp.time2.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 12-14 >> $macro/concat/$micro.hap.denovo.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 12-14 >> $macro/concat/$micro.snp.denovo.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/mathieson.txt >> $macro/concat/$micro.mathieson.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/slimulation.freq >> $macro/concat/$micro.slim.freq.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/clues.v1.selcoef >> $macro/concat/$micro.clues.txt ; done ;
# paste files
paste $macro/concat/$micro.second.mode.txt $macro/concat/$micro.first.mode.txt $macro/concat/$micro.hap.third.mode.txt $macro/concat/$micro.snp.third.mode.txt > $macro/forplots/$micro.modes.tsv
paste $macro/concat/$micro.hap.freq.txt $macro/concat/$micro.snp.freq.txt $macro/concat/$micro.slim.freq.txt > $macro/forplots/$micro.freq.tsv
paste $macro/concat/$micro.hap.sel.tsv $macro/concat/$micro.snp.sel.tsv $macro/concat/$micro.mathieson.txt $macro/concat/$micro.clues.txt | cut -f 1,4,7,10 > $macro/forplots/$micro.sel.plots.tsv
paste $macro/concat/$micro.hap.sel.tsv $macro/concat/$micro.snp.sel.tsv $macro/concat/$micro.mathieson.txt > $macro/fortables/$micro.sel.tables.tsv
paste $macro/concat/$micro.hap.time5.tsv $macro/concat/$micro.hap.time2.tsv $macro/concat/$micro.hap.denovo.tsv > $macro/fortables/$micro.hap.time.tables.tsv
paste $macro/concat/$micro.snp.time5.tsv $macro/concat/$micro.snp.time2.tsv $macro/concat/$micro.snp.denovo.tsv > $macro/fortables/$micro.snp.time.tables.tsv

micro=LARGEsLATEt
left=1
right=100
# write files
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/first.pos.txt | cut -f 2 >> $macro/concat/$micro.second.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/first.mode.txt | cut -f 2 >> $macro/concat/$micro.first.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do head -n 1 $macro/$micro/$j/third.best.hap.txt | cut -f 2 >> $macro/concat/$micro.hap.third.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do head -n 1 $macro/$micro/$j/third.best.snp.txt | cut -f 2 >> $macro/concat/$micro.snp.third.mode.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 2-4 >> $macro/concat/$micro.hap.sel.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 2-4 >> $macro/concat/$micro.snp.sel.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 1 >> $macro/concat/$micro.hap.freq.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 1 >> $macro/concat/$micro.snp.freq.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 6-8 >> $macro/concat/$micro.hap.time5.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 6-8 >> $macro/concat/$micro.snp.time5.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 9-11 >> $macro/concat/$micro.hap.time2.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 9-11 >> $macro/concat/$micro.snp.time2.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.hap.tsv | cut -f 12-14 >> $macro/concat/$micro.hap.denovo.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/results.snp.tsv | cut -f 12-14 >> $macro/concat/$micro.snp.denovo.tsv ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/mathieson.txt >> $macro/concat/$micro.mathieson.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/slimulation.freq >> $macro/concat/$micro.slim.freq.txt ; done ;
for j in $(seq $left 1 $right) ; do tail -n 1 $macro/$micro/$j/clues.v1.selcoef >> $macro/concat/$micro.clues.txt ; done ;
# paste files
paste $macro/concat/$micro.second.mode.txt $macro/concat/$micro.first.mode.txt $macro/concat/$micro.hap.third.mode.txt $macro/concat/$micro.snp.third.mode.txt > $macro/forplots/$micro.modes.tsv
paste $macro/concat/$micro.hap.freq.txt $macro/concat/$micro.snp.freq.txt $macro/concat/$micro.slim.freq.txt > $macro/forplots/$micro.freq.tsv
paste $macro/concat/$micro.hap.sel.tsv $macro/concat/$micro.snp.sel.tsv $macro/concat/$micro.mathieson.txt $macro/concat/$micro.clues.txt | cut -f 1,4,7,10 > $macro/forplots/$micro.sel.plots.tsv
paste $macro/concat/$micro.hap.sel.tsv $macro/concat/$micro.snp.sel.tsv $macro/concat/$micro.mathieson.txt > $macro/fortables/$micro.sel.tables.tsv
paste $macro/concat/$micro.hap.time5.tsv $macro/concat/$micro.hap.time2.tsv $macro/concat/$micro.hap.denovo.tsv > $macro/fortables/$micro.hap.time.tables.tsv
paste $macro/concat/$micro.snp.time5.tsv $macro/concat/$micro.snp.time2.tsv $macro/concat/$micro.snp.denovo.tsv > $macro/fortables/$micro.snp.time.tables.tsv

#!/bin/bash
#
#SBATCH --time=1-00
#SBATCH --nodes=1
#SBATCH --ntasks=1 # Number of cores
#SBATCH --mem=20000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production # Partition to submit to
#SBATCH --account=epigenetics
#SBATCH --reservation=epigenetics-workshop # reservation to submit to
#SBATCH --output=slurmout/dt-%j.out # File to which STDOUT will be written
#SBATCH --error=slurmout/dt-%j.err # File to which STDERR will be written

start=`date +%s`
hostname

aklog

export baseP=/share/workshop/epigenetics_workshop/$USER/Methylation
export inP=${baseP}/04-Methylation
export outP=${baseP}/05-Deeptools
export seqP=${baseP}/02-Cleaned
export refP=${baseP}/References
export cwd=${baseP}/scripts
export tmpd=$cwd/tmp


if [ ! -d "${outP}" ]; then
   mkdir -p ${outP}
fi

if [ ! -d "${bamP}" ]; then
   mkdir -p ${bamP}
fi

module load deeptools/3.3.1
module load kentutils/302.0.0

ID1=`head -n 1 allsamples.txt |tail -1`
ID2=`head -n 2 allsamples.txt |tail -1`
ID3=`head -n 3 allsamples.txt |tail -1`

for name in ${ID1} ${ID2} ${ID3}
do
  zcat ${inP}/${name}.input${name}_R1_001_bismark_bt2_pe.deduplicated.bedGraph.gz |tail -n +2 - |sort -k1,1 -k2,2n - > ${inP}/${name}.input${name}_R1_001_bismark_bt2_pe.deduplicated.sorted.bedGraph

  bedGraphToBigWig ${inP}/${name}.input${name}_R1_001_bismark_bt2_pe.deduplicated.sorted.bedGraph $refP/chr18.fa.fai ${outP}/${name}.input${name}_R1_001_bismark_bt2_pe.deduplicated.bw

done

computeMatrix reference-point -S ${outP}/${ID1}.input${ID1}_R1_001_bismark_bt2_pe.deduplicated.bw ${outP}/${ID2}.input${ID2}_R1_001_bismark_bt2_pe.deduplicated.bw ${outP}/${ID3}.input${ID3}_R1_001_bismark_bt2_pe.deduplicated.bw -R ../References/transcripts.noexpressed.0.bed --referencePoint TSS -a 2000 -b 2000 --outFileName ${outP}/all.matrix_Genes_noexpressed.0.tab.gz --outFileNameMatrix ${outP}/all.matrix_Genes_noexpressed.0.tab

computeMatrix reference-point -S ${outP}/${ID1}.input${ID1}_R1_001_bismark_bt2_pe.deduplicated.bw ${outP}/${ID2}.input${ID2}_R1_001_bismark_bt2_pe.deduplicated.bw ${outP}/${ID3}.input${ID3}_R1_001_bismark_bt2_pe.deduplicated.bw -R ../References/transcripts.expressed.0.bed --referencePoint TSS -a 2000 -b 2000 --outFileName ${outP}/all.matrix_Genes_expressed.0.tab.gz --outFileNameMatrix ${outP}/all.matrix_Genes_expressed.0.tab


plotHeatmap -m ${outP}/all.matrix_Genes_noexpressed.0.tab.gz -out ${outP}/all.Genes_noexpressed.0.png --colorMap jet --missingDataColor "#FFF6EB" --heatmapHeight 15 --yMin 0 --yMax 100 --plotTitle 'Not expressed genes' --regionsLabel 'genes not expressed'

plotHeatmap -m ${outP}/all.matrix_Genes_expressed.0.tab.gz -out ${outP}/all.Genes_expressed.0.png --colorMap jet --missingDataColor "#FFF6EB" --heatmapHeight 15 --yMin 0 --yMax 100 --plotTitle 'Expressed genes' --regionsLabel 'Genes expressed'


end=`date +%s`
runtime=$((end - start ))
echo $runtime



#!/bin/bash
#
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1 # Number of cores
#SBATCH --mem=20000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production # Partition to submit to
#SBATCH --account=epigenetics
#SBATCH --reservation=epigenetics-workshop
#SBATCH --output=slurmout/hm-%N-%j.out # File to which STDOUT will be written
#SBATCH --error=slurmout/hm-%N-%j.err # File to which STDERR will be written

start=`date +%s`
hostname

aklog

module load homer/4.10

export baseP=/share/workshop/epigenetics_workshop/$USER/Methylation
export outP=$baseP/05-DM
export refP=$baseP/References


tail -n +2 $outP/Differential_methylation_regions_24Mvs3M.txt |awk 'BEGIN{FS="\t";OFS="\t"}{print $1, $2 - 1, $3}' - |sort -k1,1 -k2,2n - > $outP/Differential_methylation_regions_24Mvs3M.bed


annotatePeaks.pl $outP/Differential_methylation_regions_24Mvs3M.bed mm10 -raw -gid -annStats $outP/homer.stats -go $outP/GO > $outP/Differential_methylation_regions_24Mvs3M.anno.txt


end=`date +%s`
runtime=$((end - start ))
echo $runtime



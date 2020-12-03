#!/bin/bash
#
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1 # Number of cores
#SBATCH --mem=10000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production # Partition to submit to
#SBATCH --account=epigenetics
#SBATCH --reservation=epigenetics-workshop
#SBATCH --output=slurmout/bt-%N-%j.out # File to which STDOUT will be written
#SBATCH --error=slurmout/bt-%N-%j.err # File to which STDERR will be written

start=`date +%s`
hostname

aklog

module load bedtools2/2.29.2

export baseP=/share/workshop/epigenetics_workshop/$USER/Methylation
export outP=$baseP/05-DM
export refP=$baseP/References




bedtools intersect -wo -a $outP/Differential_methylation_regions_24Mvs3M.bed -b $refP/GSM1865011_ATAC-seq_WT_rods_MACS_peaks_rep1.txt > $outP/DMR.ATAC.txt


end=`date +%s`
runtime=$((end - start ))
echo $runtime



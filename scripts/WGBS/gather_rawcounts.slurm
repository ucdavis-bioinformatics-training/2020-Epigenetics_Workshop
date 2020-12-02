#!/bin/bash
#
#
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1 # Number of cores
#SBATCH --mem=10g # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production # Partition to submit to
#SBATCH --output=slurmout/cnt-%A-%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/cnt-%A-%a.err # File to which STDERR will be written

start=`date +%s`
hostname

aklog


export baseP=/share/workshop/epigenetics_workshop/$USER/Methylation
export bamP=${baseP}/03-Bismark
export outP=${baseP}/04-Methylation
export seqP=${baseP}/02-Cleaned
export refP=${baseP}/References
export sxpd=${baseP}/scripts


if [ ! -d "${outP}" ]; then
   mkdir -p ${outP}
fi

if [ ! -d "${bamP}" ]; then
   mkdir -p ${bamP}
fi


name=`head -n ${SLURM_ARRAY_TASK_ID} ${sxpd}/allsamples.txt | tail -1`
input=`ls $outP/${name}.input*.bismark.cov.gz`

# bismark bedGraph output coverage file: <chromosome> <start> <end> <methylation percentage> <count methylated> <count unmethylated>
# DSS methylation count input file: <chromosome> <pos> <count total> <count methylated>
zcat ${input} |awk 'BEGIN{FS="\t";OFS="\t"}{print $1, $2, ($5 + $6), $5}' - > $outP/${name}.CpG.counts.txt


end=`date +%s`
runtime=$((end - start ))
echo $runtime



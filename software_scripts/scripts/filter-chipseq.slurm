#!/bin/bash
#
#SBATCH --job-name=filter-chip # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=1 # Number of cores
#SBATCH --mem=4000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --time=6:00:00
#SBATCH --array=1-8
#SBATCH --partition=production # Partition to submit to
#SBATCH --account=epigenetics # cluster account to use for the job
#SBATCH --reservation=epigenetics-workshop # cluster account reservation
#SBATCH --output=slurm_out/filter-%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurm_out/filter-%A_%a.err # File to which STDERR will be written

start=`date +%s`
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

inpath=02-BWA
sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" samples.txt | awk -F '\t'  '{print $1}'`
bam=${inpath}/${sample}/${sample}_bwa.bam

outpath='03-Filter'
[[ -d ${outpath} ]] || mkdir ${outpath}
[[ -d ${outpath}/${sample} ]] || mkdir ${outpath}/${sample}

echo "SAMPLE: ${sample}"

module load samtools/1.9
module load deeptools/3.3.1
module load bedtools2/2.29.2

output=${outpath}/${sample}/${sample}

THREADS=${SLURM_NTASKS}
#THREADS=1

blacklist_bed=/share/workshop/epigenetics_workshop/msettles/References/mm10-blacklist.v2.bed

call="samtools view -F 2048 -f 2 -q 30 -b ${bam} | samtools sort -m 768M --threads ${THREADS} -o ${output}_filtered.bam -"
echo $call
eval $call

call="samtools index -@ ${THREADS} ${output}_filtered.bam"
echo $call
eval $call

call="bedtools intersect -abam ${output}_filtered.bam -b ${blacklist_bed} -v | samtools sort -m 768M --threads ${THREADS} -o ${output}_filtered_blacklisted.bam -"
echo $call
eval $call

call="samtools index -@ ${THREADS} ${output}_filtered_blacklisted.bam"
echo $call
eval $call

call="samtools idxstats ${output}_filtered_blacklisted.bam > ${output}_filtered_blacklisted.bam.idxstats"
echo $call
eval $call

call="samtools flagstat -@ ${THREADS} ${output}_filtered_blacklisted.bam > ${output}_filtered_blacklisted.bam.flagstat"
echo $call
eval $call

call="samtools stats -@ ${THREADS} ${output}_filtered_blacklisted.bam > ${output}_filtered_blacklisted.bam.stats"
echo $call
eval $call

end=`date +%s`

runtime=$((end-start))

echo $runtime

#!/bin/bash
#
#SBATCH --job-name=map_hts_bwa # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=20 # Number of cores
#SBATCH --mem=16000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --time=1-00
#SBATCH --array=1-12
#SBATCH --partition=production # Partition to submit to
#SBATCH --account=epigenetics # cluster account to use for the job
#SBATCH --reservation=epigenetics-workshop # cluster account reservation
#SBATCH --output=slurm_out/map_hts_bwa-%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurm_out/map_hts_bwa-%A_%a.err # File to which STDERR will be written
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=settles@ucdavis.edu # Email to which notifications will be sent

start=`date +%s`
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

inpath=01-HTS_Preproc
sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" samples.txt | awk -F '\t'  '{print $1}'`
r1=${inpath}/${sample}/${sample}*_R1*.fastq.gz
r2=${inpath}/${sample}/${sample}*_R2*.fastq.gz

outpath='02-BWA'
[[ -d ${outpath} ]] || mkdir ${outpath}
[[ -d ${outpath}/${sample} ]] || mkdir ${outpath}/${sample}

echo "SAMPLE: ${sample}"

THREADS=${SLURM_NTASKS}
MEM=$(expr ${SLURM_MEM_PER_NODE} / 1024)
MAPTHREADS=$(expr ${THREADS} - 6)
SORTTHREADS=$(expr ${THREADS} - ${MAPTHREADS})

#THREADS=20
#MAPTHREADS=14
#SORTTHREADS=6
#MEM=48

module load bwa/0.7.17
module load samtools/1.9

output=${outpath}/${sample}/${sample}_bwa.bam
mapfasta=/share/biocore/workshops/2020_Epigenetics/Reference/GRCm38.primary_assembly.genome.fa

call="bwa mem -t ${MAPTHREADS} \
 -R '@RG\tID:${sample}\tSM:${new_id}\tPL:ILLUMINA\tDS:Paired' \
 ${mapfasta} ${r1} ${r2} \
 | samtools sort -m 768M --threads ${SORTTHREADS} -o ${output} -"
echo $call
eval $call

call="samtools index -@ ${THREADS} ${output}"
echo $call
eval $call

call="samtools idxstats ${output} > ${output}.idxstats"
echo $call
eval $call

call="samtools flagstat -@ ${THREADS} ${output} > ${output}.flagstat"
echo $call
eval $call

call="samtools stats -@ ${THREADS} ${output} > ${output}.stats"
echo $call
eval $call

end=`date +%s`

runtime=$((end-start))

echo $runtime

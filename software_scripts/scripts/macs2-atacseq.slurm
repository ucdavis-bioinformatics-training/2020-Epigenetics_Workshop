#!/bin/bash
#
#SBATCH --job-name=macs-atac # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=1 # Number of cores
#SBATCH --mem=4000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --time=2:00:00
#SBATCH --array=1-6
#SBATCH --partition=production # Partition to submit to
#SBATCH --account=epigenetics # cluster account to use for the job
#SBATCH --reservation=epigenetics-workshop # cluster account reservation
#SBATCH --output=slurm_out/macs2-%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurm_out/macs2-%A_%a.err # File to which STDERR will be written
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=settles@ucdavis.edu # Email to which notifications will be sent

start=`date +%s`
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

inpath=03-Filter
sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" samples.txt | awk -F '\t'  '{print $1}'`
bam=${inpath}/${sample}/${sample}_shifted_filtered_blacklisted.bam

outpath='04-MACS2'
[[ -d ${outpath} ]] || mkdir ${outpath}
[[ -d ${outpath}/${sample} ]] || mkdir ${outpath}/${sample}

echo "SAMPLE: ${sample}"

output=${outpath}/${sample}/${sample}

module load  macs2/2.2.5

THREADS=${SLURM_NTASKS}
#THREADS=1

call="macs2 callpeak  -t ${bam}  -f BAMPE  -n ${output}  -g mm  --keep-dup all 2> ${output}.log"
echo $call
eval $call

end=`date +%s`

runtime=$((end-start))

echo $runtime

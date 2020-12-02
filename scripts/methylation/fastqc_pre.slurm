#!/bin/bash
#
#SBATCH --time=3-00
#SBATCH --nodes=1
#SBATCH --ntasks=8 # Number of cores
#SBATCH --mem=80000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production # Partition to submit to
#SBATCH --account=epigenetics
#SBATCH --reservation=epigenetics-workshop # reservation to submit to
#SBATCH --output=slurmout/fqc-%A-%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/fqc-%A-%a.err # File to which STDERR will be written

start=`date +%s`
hostname

aklog
export baseP=/share/workshop/epigenetics_workshop/$USER/Methylation
export seqP=${baseP}/00-RawData
export outP=${baseP}/01-Fastqc
export cwd=${baseP}/scripts
export tmpP=$cwd/tmp


if [ ! -d "${outP}" ]; then
   mkdir ${outP}
fi

if [ ! -d "${tmpP}" ]; then
   mkdir ${tmpP}
fi


module load fastqc/0.11.9

name=`head -n ${SLURM_ARRAY_TASK_ID} allsamples.txt | tail -1`


call="fastqc -t ${SLURM_NTASKS} --dir $tmpP --outdir $outP $seqP/${name}_1.fastq"
echo $call
eval $call

call="fastqc -t ${SLURM_NTASKS} --dir $tmpP --outdir $outP $seqP/${name}_2.fastq"
echo $call
eval $call


end=`date +%s`
runtime=$((end - start ))
echo $runtime



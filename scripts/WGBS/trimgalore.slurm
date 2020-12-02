#!/bin/bash
#
#SBATCH --time=1-00
#SBATCH --nodes=1
#SBATCH --ntasks=24 # Number of cores
#SBATCH --mem=80000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production # Partition to submit to
#SBATCH --output=slurmout/tg-%A-%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/tg-%A-%a.err # File to which STDERR will be written

start=`date +%s`
hostname

aklog
export baseP=/share/workshop/epigenetics_workshop/$USER/Methylation
export seqP=${baseP}/00-RawData
export outP=${baseP}/02-Cleaned
export cwd=${baseP}/scripts


if [ ! -d "${outP}" ]; then
   mkdir ${outP}
fi


module load samtools/1.10
module load bowtie2/2.4.2
module load trim_galore/0.6.6
module load cutadapt/2.10

source activate cutadapt-2.10

name=`head -n ${SLURM_ARRAY_TASK_ID} allsamples.txt | tail -1`


  R1=$(ls $seqP/${name}_*.fastq |grep "_1." - )
  R2=$(ls $seqP/${name}_*.fastq |grep "_2." - )

  call="trim_galore -j ${SLURM_NTASKS} --paired -q 30 --stringency 5 --phred33 -e 0.1 --length 20 \
        --clip_R1 10 --clip_R2 20 --three_prime_clip_R1 10 --three_prime_clip_R2 15 \
        --output_dir $outP --basename ${name} ${R1} ${R2}"
  echo $call
  eval $call


end=`date +%s`
runtime=$((end - start ))
echo $runtime



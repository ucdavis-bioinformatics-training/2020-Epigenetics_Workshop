#!/bin/bash
#
#SBATCH --time=1-00
#SBATCH --nodes=1
#SBATCH --ntasks=1 # Number of cores
#SBATCH --mem=10000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production # Partition to submit to
#SBATCH --output=slurmout/gp-%N-%j.out # File to which STDOUT will be written
#SBATCH --error=slurmout/gp-%N-%j.err # File to which STDERR will be written

export baseP=/share/workshop/epigenetics_workshop/$USER/Methylation
export bamP=${baseP}/03-Bismark
export seqP=${baseP}/02-Cleaned
export refP=${baseP}/References
export cwd=${baseP}/scripts


if [ ! -d "${outP}" ]; then
   mkdir ${outP}
fi

if [ ! -d "${bamP}" ]; then
   mkdir ${bamP}
fi

module load bismark/0.22.3
module load bowtie2/2.4.2

bismark_genome_preparation --bowtie2 $refP


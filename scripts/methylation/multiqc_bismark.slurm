#!/bin/bash
#
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1 # Number of cores
#SBATCH --mem=10000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production # Partition to submit to
#SBATCH --account=epigenetics
#SBATCH --reservation=epigenetics-workshop # reservation to submit to
#SBATCH --output=slurmout/mqcb-%j.out # File to which STDOUT will be written
#SBATCH --error=slurmout/mqcb-%j.err # File to which STDERR will be written

start=`date +%s`
hostname

aklog
export baseP=/share/workshop/epigenetics_workshop/$USER/Methylation
export alnP=${baseP}/03-Bismark
export metP=${baseP}/04-Methylation
export cwd=${baseP}/scripts
export tmpP=$cwd/tmp


if [ ! -d "${alnP}" ]; then
   mkdir -p ${alnP}
fi

if [ ! -d "${tmpP}" ]; then
   mkdir -p ${tmpP}
fi


module load multiqc/1.9


call="multiqc -l input2.fofn -m bismark -o $metP/"
echo $call
eval $call


end=`date +%s`
runtime=$((end - start ))
echo $runtime



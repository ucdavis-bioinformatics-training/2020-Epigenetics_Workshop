#!/bin/bash
#
#SBATCH --time=1-00
#SBATCH --nodes=1
#SBATCH --ntasks=2 # Number of cores
#SBATCH --mem=20000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production # Partition to submit to
#SBATCH --output=slurmout/mqc-%j.out # File to which STDOUT will be written
#SBATCH --error=slurmout/mqc-%j.err # File to which STDERR will be written

start=`date +%s`
hostname

aklog
export baseP=/share/workshop/epigenetics_workshop/$USER/Methylation
export seqP=${baseP}/00-RawData
export outP=${baseP}/011-Fastqc
export cwd=${baseP}/scripts
export tmpP=$cwd/tmp


if [ ! -d "${outP}" ]; then
   mkdir ${outP}
fi

if [ ! -d "${tmpP}" ]; then
   mkdir ${tmpP}
fi


module load multiqc/1.9


call="multiqc -f -o $outP/ $outP/*_fastqc.zip"
echo $call
eval $call


end=`date +%s`
runtime=$((end - start ))
echo $runtime



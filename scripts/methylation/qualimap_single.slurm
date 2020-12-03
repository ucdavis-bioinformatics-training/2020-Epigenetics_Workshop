#!/bin/bash
#SBATCH --time=1-00
#SBATCH --nodes=1
#SBATCH --ntasks=4 # Number of cores
#SBATCH --mem=50000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production # Partition to submit to
#SBATCH --account=epigenetics
#SBATCH --reservation=epigenetics-workshop # reservation to submit to
#SBATCH --output=slurmout/qm-%A-%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/qm-%A-%a.err # File to which STDERR will be written


start=`date +%s`
hostname

export baseP=/share/workshop/epigenetics_workshop/$USER/Methylation
export bamP=${baseP}/03-Bismark
export seqP=${baseP}/02-Cleaned
export refP=${baseP}/References
export cwd=${baseP}/scripts


if [ ! -d "${bamP}" ]; then
   mkdir ${bamP}
fi

module load samtools/1.10
module load qualimap/2.2.1

name=`head -n ${SLURM_ARRAY_TASK_ID} allsamples.txt | tail -1`

export outP=${baseP}/04-Qualimap/${name}

if [ ! -d "${outP}" ]; then
   mkdir -p ${outP}
fi

samtools sort -@ ${SLURM_NTASKS} -m 10G -o $bamP/${name}.input${name}_R1_001_bismark_bt2_pe.dedup.sorted.bam $bamP/${name}.input${name}_R1_001_bismark_bt2_pe.deduplicated.bam

qualimap bamqc --java-mem-size=40G -bam ${bamP}/${name}.input${name}_R1_001_bismark_bt2_pe.dedup.sorted.bam -nt ${SLURM_NTASKS} -c -oc $outP/${name}.qualimap.txt -outdir ${outP} -outfile ${name}.html -outformat html



end=`date +%s`
runtime=$((end - start ))
echo $runtime



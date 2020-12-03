#!/bin/bash
#SBATCH --time=1:00
#SBATCH --nodes=1
#SBATCH --ntasks=2 # Number of cores
#SBATCH --mem=10000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production # Partition to submit to
#SBATCH --account=epigenetics
#SBATCH --reservation=epigenetics-workshop # reservation to submit to
#SBATCH --output=slurmout/qmm-%j.out # File to which STDOUT will be written
#SBATCH --error=slurmout/qmm-%j.err # File to which STDERR will be written


start=`date +%s`
hostname

export baseP=/share/workshop/epigenetics_workshop/$USER/Methylation
export datP=${baseP}/04-Qualimap
export seqP=${baseP}/02-Cleaned
export refP=${baseP}/References
export cwd=${baseP}/scripts


module load samtools/1.10
module load qualimap/2.2.1

n=$(wc -l allsamples.txt |cut -d' ' -f1 )

for i in $( eval echo {1..${n}} )
do
  sample=$(sed "${i}q;d" allsamples.txt)
  if [ $i -lt 4 ]
  then
    grp="group_1"
  else
    grp="group_2"
  fi
  echo -e "${sample}\t${datP}/${sample}\t${grp}" >> input.fofn
done


qualimap multi-bamqc -d input.fofn -outdir $datP 


end=`date +%s`
runtime=$((end - start ))
echo $runtime



#!/bin/bash

#SBATCH --job-name=htstream # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=9
#SBATCH --time=12:00:00
#SBATCH --mem=15000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --array=1-8
#SBATCH --partition=production
#SBATCH --account=epigenetics # cluster account to use for the job
#SBATCH --reservation=epigenetics-workshop # cluster account reservation
#SBATCH --output=slurm_out/htstream_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurm_out/htstream_%A_%a.err # File to which STDERR will be written

start=`date +%s`
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

inpath="00-RawData"

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" samples.txt | awk -F '\t'  '{print $1}'`
r1=${inpath}/${sample}/${sample}*_R1*.fastq.gz
r2=${inpath}/${sample}/${sample}*_R2*.fastq.gz

outpath='01-HTS_Preproc'
[[ -d ${outpath} ]] || mkdir ${outpath}
[[ -d ${outpath}/${sample} ]] || mkdir ${outpath}/${sample}

echo "SAMPLE: ${sample}"

module load htstream/1.3.2

call="hts_Stats -L ${outpath}/${sample}/${sample}_htsStats.log -1 ${r1} -2 ${r2} -N 'initial Stats' | \
      hts_SeqScreener -A  ${outpath}/${sample}/${sample}_htsStats.log -N 'PhiX check' | \
      hts_SuperDeduper -e 250000 -A  ${outpath}/${sample}/${sample}_htsStats.log -N 'Remove PCR duplicates' | \
      hts_AdapterTrimmer -p 4 -A  ${outpath}/${sample}/${sample}_htsStats.log -N 'Overlap and remove adapters' | \
      hts_NTrimmer -A ${outpath}/${sample}/${sample}_htsStats.log -N 'Remove all Ns' | \
      hts_QWindowTrim -A ${outpath}/${sample}/${sample}_htsStats.log -N 'Quality trim' | \
      hts_LengthFilter -n -m 50 -A ${outpath}/${sample}/${sample}_htsStats.log -N 'Remove too short' | \
      hts_Stats -A ${outpath}/${sample}/${sample}_htsStats.log -F -f ${outpath}/${sample}/${sample} -N 'end Stats'"

echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime

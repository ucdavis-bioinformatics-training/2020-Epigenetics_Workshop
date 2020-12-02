#!/bin/bash

#SBATCH --job-name=bam_index # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=120
#SBATCH --mem=40000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production
#SBATCH --account=epigenetics # cluster account to use for the job
#SBATCH --reservation=epigenetics-workshop # cluster account reservation
#SBATCH --output=bwa-index_%A.out # File to which STDOUT will be written
#SBATCH --error=bwa-index_%A.err # File to which STDERR will be written

start=`date +%s`
echo $HOSTNAME

outpath="References"
mkdir -p ${outpath}

cd ${outpath}

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz GRCm38.primary_assembly.genome.fa.gz
gunzip GRCm38.primary_assembly.genome.fa.gz
FASTA="GRCm38.primary_assembly.genome.fa"

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.primary_assembly.annotation.gtf.gz gencode.vM25.primary_assembly.annotation.gtf.gz
gunzip gencode.vM25.primary_assembly.annotation.gtf.gz


module load bwa/0.7.17

call="bwa index ${FASTA}"

echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime

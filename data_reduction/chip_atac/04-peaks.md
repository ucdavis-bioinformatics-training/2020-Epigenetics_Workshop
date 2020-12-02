
# Calling Peaks using MACS2

This document assumes [filtering of samples](./03-filter.md) has been completed.

**IF** for some reason it didn't finish, is corrupted or you missed the session, you can copy over a completed copy

```bash
#cp -r /share/biocore/workshops/2020_Epigenetics/ChIPseq/03-Filter /share/workshop/epigenetics_workshop/$USER/chipseq_example/.
#cp -r /share/biocore/workshops/2020_Epigenetics/ATACseq/03-Filter /share/workshop/epigenetics_workshop/$USER/atacseq_example/.
```

# Peak Calling with MACS2 for ChIPseq


1. We can now run MACS2 Peak calling across all samples on the real data using a SLURM script, [macs2-chipseq.slurm](../../software_scripts/scripts/macs2-chipseq.slurm), that we should take a look at now.

    ```bash
    cd /share/workshop/epigenetics_workshop/$USER/chipseq_example  # We'll run this from the main directory

    wget https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2020-Epigenetics_Workshop/master/software_scripts/scripts/macs2-chipseq.slurm macs2-chipseq.slurm
    less macs2-chipseq.slurm
    ```

    <div class="script">#!/bin/bash
    #
    #SBATCH --job-name=macs_chip # Job name
    #SBATCH --nodes=1
    #SBATCH --ntasks=1 # Number of cores
    #SBATCH --mem=4000 # Memory pool for all cores (see also --mem-per-cpu)
    #SBATCH --time=2:00:00
    #SBATCH --array=1-8
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
    bam=${inpath}/${sample}/${sample}_filtered_blacklisted.bam
    input="JLDY037L"
    input_bam=${inpath}/${input}/${input}_filtered_blacklisted.bam

    outpath='04-MACS2'
    [[ -d ${outpath} ]] || mkdir ${outpath}
    [[ -d ${outpath}/${sample} ]] || mkdir ${outpath}/${sample}

    echo "SAMPLE: ${sample}"

    output=${outpath}/${sample}/${sample}

    module load  macs2/2.2.5

    THREADS=${SLURM_NTASKS}
    #THREADS=1

    #call="macs2 callpeak  -t ${bam}  -c ${input_bam} -f BAMPE  -n ${output}  -g mm  --keep-dup all 2> ${output}.log"
    call="macs2 callpeak  -t ${bam} -f BAMPE  -n ${output}  -g mm  --keep-dup all 2> ${output}.log"
    echo $call
    eval $call

    end=`date +%s`

    runtime=$((end-start))

    echo $runtime
    </div>

2. After looking at the script, lets run it.

    ```bash
    sbatch macs2-atacseq.slurm  # moment of truth!
    ```

    We can watch the progress of our task array using the 'squeue' command. Takes about 1:30 hours to process each sample.

    ```sbatch
    squeue -u $USER  # use your username
    ```

## MACS2 for ATAC

1. We can now run MACS2 Peak calling across all samples on the real data using a SLURM script, [macs2-atacseq.slurm](../../software_scripts/scripts/macs2-atacseq.slurm), that we should take a look at now.

    ```bash
    cd /share/workshop/epigenetics_workshop/$USER/atacseq_example  # We'll run this from the main directory

    wget https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2020-Epigenetics_Workshop/master/software_scripts/scripts/macs2-atacseq.slurm macs2-atac.slurm
    less macs2-atacseq.slurm
    ```

    <div class="script">#!/bin/bash
    #
    #SBATCH --job-name=macs-atac # Job name
    #SBATCH --nodes=1
    #SBATCH --ntasks=1 # Number of cores
    #SBATCH --mem=4000 # Memory pool for all cores (see also --mem-per-cpu)
    #SBATCH --time=2
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
    </div>

2. After looking at the script, lets run it.

    ```bash
    sbatch macs2-atacseq.slurm  # moment of truth!
    ```

    We can watch the progress of our task array using the 'squeue' command. Takes about 1:30 hours to process each sample.

    ```sbatch
    squeue -u $USER  # use your username
    ```

## MultiQC QA/QC Summary of the filter results.

Well there is a [MultiQC](https://multiqc.info/) report for MACS but its pretty pointless to generate.

```bash
## Run multiqc to collect statistics and create a report:
cd /share/workshop/epigenetics_workshop/$USER/chipseq_example
module load multiqc/htstream.dev0
multiqc -i ChIPseq-macs2-report -o 04-MACS2-ChIPseq-report ./04-MACS2
```

**Do the same for the ATACseq experiment**


Transfer ChIPseq-macs-report_multiqc_report.html and ATACseq-filter-report_multiqc_report.html to your computer and open it in a web browser.

Or in case of emergency, download this copy: [ChIPseq-macs2-report_multiqc_report.html](ChIPseq-macs2-report_multiqc_report.html) and [ATACseq-macs2-report_multiqc_report.html](ATACseq-macs2-report_multiqc_report.html) for the ATACseq

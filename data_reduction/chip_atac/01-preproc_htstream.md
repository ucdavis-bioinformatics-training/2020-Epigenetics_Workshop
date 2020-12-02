# Sequence Preprocessing

This document assumes [project_setup](./00-project_setup.md) has been completed.

```bash
cd /share/workshop/epigenetics_workshop/$USER/chipseq_example
```

## Why Preprocess Reads

We have found that aggressively “cleaning” and preprocessing of reads can make a large difference to the speed and quality of mapping and assembly results. Preprocessing your reads means:

  * Removing bases of unwanted sequence (Ex. vectors, adapter, primer sequence, polyA tails).
  * Merge/join short overlapping paired-end reads.
  * Remove low quality bases or N characters.
  * Remove reads originating from PCR duplication.
  * Remove reads that are not of primary interest (contamination).
  * Remove too short reads.

Preprocessing also produces a number of statistics about the samples. These can be used to evaluate sample-to-sample consistency.

### Preprocessing Statistics as QA/QC.

Beyond generating "better" data for downstream analysis, preprocessing statistics also give you an idea as to the original quality and complexity of the sample, library generation features, and sequencing quality.

This can help inform you of how you might change your procedures in the future, either sample preparation, or in library preparation.

We’ve found it best to perform __QA/QC__ on both the run as a whole (poor samples can negatively affect other samples) and on the samples themselves as they compare to other samples (**BE CONSISTENT**).

Reports such as Basespace for Illumina, are great ways to evaluate the run as a whole, the sequencing provider usually does this for you.  

Visualizing the preprocessing summary are a great way to look for technical bias across your experiment. Poor quality samples often appear as outliers and can ethically be removed due to identified technical issues. You should **NOT** see a trend associated with any experimental factors. That scenario should raise concerns about technical sample processing bias.

**Many technical things happen between original sample and data, preprocessing is working backwards through that process to get as close as we can to original sample**

<img src="preproc_figures/preproc_flowchart.png" alt="preproc_flowchart" width="80%"/>

### An ChIPseq/ATACseq Preprocessing Workflow

1. Raw data stats.
1. Remove contaminants (at least PhiX).
1. Remove PCR duplicates.
1. Overlapping paired end reads and remove adapters (overhangs).
1. Trim off all 'N' bases.
1. Trim sequences (5’ and 3’) by quality score (I like Q20).
1. Cleanup.
  * Remove any reads that are less then the minimum length parameter.
  * Produce preprocessing statistics.

## HTStream Streamed Preprocessing of Sequence Data

HTStream is a suite of preprocessing applications for high throughput sequencing data (ex. Illumina). A fast C++ implementation, designed with discreet functionality that can be pipelined together using standard Unix piping.

Benefits Include:
  * No intermediate files, reducing storage footprint.
  * Reduced I/O, files are only read in and written out once to disk.
  * Handles both single end and paired end reads at the same time.
  * Applications process reads at the same time allowing for process parallelization.
  * Built on top of mature C++ Boost libraries to reduce bugs and memory leaks.
  * Designed following the philosophy of [Program Design in the UNIX Environment](https://onlinelibrary.wiley.com/doi/abs/10.1002/j.1538-7305.1984.tb00055.x).
  * Works with native Unix/Linux applications such as grep/sed/awk etc.
  * Can build a custom preprocessing pipeline to fit the specific expectation of the data.
  * A single JSON output per sample detailing the preprocessing statistics from each application.

HTStream achieves these benefits by using a tab delimited intermediate format that allows for streaming from application to application. This streaming creates some awesome efficiencies when preprocessing HTS data and makes it fully interoperable with other standard Linux tools.

### HTStream applications

HTStream includes the following applications:

hts_AdapterTrimmer: Identify and remove adapter sequences.  
hts_CutTrim: Discreet 5' and/or 3' basepair trimming.  
hts_LengthFilter: Remove reads outside of min and/or max length.  
hts_NTrimmer: Extract the longest subsequence with no Ns.    
hts_Overlapper: Overlap paired end reads, removing adapters when present.  
hts_PolyATTrim: Identify and remove polyA/T sequence.  
hts_Primers: Identify and optionally remove 5' and/or 3' primer sequence.  
hts_QWindowTrim: 5' and/or 3' quality score base trimming using windows.  
hts_SeqScreener: Identify and remove/keep/count contaminants (default phiX).  
hts_Stats: Compute read stats.  
hts_SuperDeduper: Identify and remove PCR duplicates.  

The source code and pre-compiled binaries for Linux can be downloaded and installed [from the GitHub repository](https://github.com/s4hts/HTStream).

HTStream is also avaiable on [Bioconda](https://bioconda.github.io/), and there is even an image on [Docker Hub](https://hub.docker.com/r/dzs74/htstream).

HTStream was designed to be extensible. We continue to add new preprocessing routines and welcome contributions from collaborators.

If you encounter any bugs or have suggestions for improvement, please post them to [issues](https://github.com/s4hts/HTStream/issues).

# HTStream Setup for our Project

## Example, running HTStream

Let's run the first step of our HTStream preprocessing pipeline, which is always to gather basic stats on the read files. For now, we're only going to run one sample through the pipeline.

When building a new pipeline, it is almost always a good idea to use a small subset of the data in order to speed up development. A small sample of reads will take seconds to process and help you identify problems that may have only been apparent after hours of waiting for the full data set to process.

1. Let's start by first taking a small subsample of reads, so that our trial run through the pipeline goes really quickly.

    ```bash
    cd /share/workshop/epigenetics_workshop/$USER/chipseq_example
    mkdir HTS_testing
    cd HTS_testing
    pwd
    ```

    * *Why run ```pwd``` here?*

    Then create a small dataset.

    ```bash
    zcat ../00-RawData/JLDY037E/JLDY037E_S5_L005_R1_001.fastq.gz | head -400000 | gzip > JLDY037E.subset_R1.fastq.gz
    zcat ../00-RawData/JLDY037E/JLDY037E_S5_L005_R2_001.fastq.gz | head -400000 | gzip > JLDY037E.subset_R2.fastq.gz
    ls
    ```

    So we ```zcat``` (uncompress and send to stdout), pipe ```|```  to ```head``` (param -400000) then pipe to ```gzip``` to recompress and name our files subset.

    * *How many reads are we going to analyze in our subset?*

1. Now we'll run our first preprocessing step ```hts_Stats```, first loading the module and then looking at help.

    ```bash
    cd /share/workshop/epigenetics_workshop/$USER/chipseq_example/HTS_testing
    module load htstream
    hts_Stats --help
    ```

    * *What version of hts_Stats is loaded?*


1. Now lets run ```hts_Stats``` and look at the output.

    ```bash
    hts_Stats -1 JLDY037E.subset_R1.fastq.gz \
              -2 JLDY037E.subset_R2.fastq.gz \
              -L JLDY037E.stats.json > out.tab
    ```

    * *What happens if you run hts_Stats without piping output to out.tab?*

    * *Can you think of a way to view the output from hts_Stats in less without creating out.tab?*

    By default, all HTS apps output tab formatted files to the stdout.

    Take a look at the output (remember ```q``` quits):
    ```bash
    less out.tab
    ```

    The output was difficult to understand, lets try without line wrapping (note that you can also type ```-S``` from within ```less``` if you forget). Scroll with the arrow keys, left, right, up, and down.
    ```bash
    less -S out.tab
    ```

    And delete out.tab since we are done with it:
    ```bash
    rm out.tab
    ```

    Remember how this output looks, we will revisit it later.

1. Now lets change the command slightly.
    ```bash
    hts_Stats -1 JLDY037E.subset_R1.fastq.gz \
              -2 JLDY037E.subset_R2.fastq.gz \
              -L JLDY037E.stats.json -f JLDY037E.stats
    ```

    * *What parameters did we use, what do they do?*

    Lets take a look at the output of stats

    ```bash
    ls -lah
    ```

    <div class="output">msettles@tadpole:/share/workshop/epigenetics_workshop/msettles/chipseq_example/HTS_testing$     ls -lah
    total 32M
    drwxrwsr-x 2 msettles epigenetics    7 Nov 29 21:22 .
    drwxrwsr-x 7 msettles epigenetics    8 Nov 29 21:16 ..
    -rw-rw-r-- 1 msettles epigenetics  60K Nov 29 21:22 JLDY037E.stats.json
    -rw-rw-r-- 1 msettles epigenetics 7.2M Nov 29 21:22 JLDY037E.stats_R1.fastq.gz
    -rw-rw-r-- 1 msettles epigenetics 8.8M Nov 29 21:22 JLDY037E.stats_R2.fastq.gz
    -rw-rw-r-- 1 msettles epigenetics 7.2M Nov 29 21:20 JLDY037E.subset_R1.fastq.gz
    -rw-rw-r-- 1 msettles epigenetics 8.8M Nov 29 21:20 JLDY037E.subset_R2.fastq.gz
    </div>

    * *Which files were generated from hts\_Stats?*

1. Lets look at the file JLDY037E.stats\.json*

    ```bash
    cat JLDY037E.stats.json
    ```

    The logs generated by htstream are in [JSON](https://en.wikipedia.org/wiki/JSON) format, like a database format but meant to be readable.


## Next lets screen out PhiX, the Illumina control

1. First, view the help documentation for hts_SeqScreener

    ```bash
    cd /share/workshop/epigenetics_workshop/$USER/chipseq_example/HTS_testing
    hts_SeqScreener -h
    ```

    * *What parameters are needed to:
        1. provide a reference to hts_SeqScreener and
        2. count, and not screen occurrences?*

1. Run HTStream on the small test set.

    ```bash
    hts_SeqScreener -1 JLDY037E.subset_R1.fastq.gz \
                    -2 JLDY037E.subset_R2.fastq.gz \
                    -r -L JLDY037E.phix.json -f JLDY037E.phix
    ```

    * *Which files were generated from hts\_SeqScreener?*

    * *Lets look at the file JLDY037E.phix.json?*

    * *What do you notice about the JLDY037E.phix.json?*

    * *How many reads were identified as phix?*

### Stream multiple applications together.

The power of HTStream is the ability to stream reads through multiple programs using pipes. By streaming reads through programs, processing will be much quicker because each read is read in only once and written out only once. This approach also uses significantly less storage as there are no intermediate files. HTStream can do this by streaming a tab-delimited format called tab6.

Single end reads are 3 columns:

`read1id  read1seq  read1qual`

Paired end reads are 6 columns:

`read1id  read1seq  read1qual  read2id  read2seq  read2qual`

1. So lets first run hts_Stats and then hts_SeqScreener in a streamed fashion.

    ```bash
    cd /share/workshop/epigenetics_workshop/$USER/chipseq_example/HTS_testing

    hts_Stats -1 JLDY037E.subset_R1.fastq.gz \
              -2 JLDY037E.subset_R2.fastq.gz \
              -L JLDY037E.streamed.json |
    hts_SeqScreener -A JLDY037E.streamed.json \
              -f JLDY037E.streamed
    ```

    Note the pipe, ```|```, between the two applications!

    **Questions**
    * *What new parameters did we use here?*

    * *What parameter is SeqScreener using that specifies how reads are input?*

    * *Lets look at the file JLDY037E.streamed.json?*


## A ChIPseq preprocessing pipeline

1. hts_Stats: get stats on *input* raw reads
1. hts_SeqScreener: screen out (remove) phiX
1. hts_SuperDeduper: identify and remove PCR duplicates
1. hts_AdapterTrimmer: identify and remove adapter sequence
1. hts_NTrimmer: trim to remove any remaining N characters
1. hts_QWindowTrim: remove poor quality bases
1. hts_LengthFilter: use to remove all reads < 50bp
1. hts_Stats: get stats on *output* cleaned reads


### Why screen for phiX?

PhiX is a common control in Illumina runs, and facilities may not tell you if/when PhiX has been spiked in. Since it does not have a barcode, in theory should not be in your data.

However:
* When we know PhiX has been spiked in, we find sequence every time.
    * [update] When dual matched barcodes are used, then almost zero phiX reads are identified.
* When I know PhiX has not been spiked in, I do not find sequence

For RNAseq and variant analysis (any mapping based technique) it is not critical to remove, but for sequence assembly it is. Unless you are sequencing PhiX, it is noise, so its better safe than sorry to screen for it every time.

### Removing PCR duplicates with hts_SuperDeduper.

Removing PCR duplicates can be **controversial** for RNAseq, but I'm in favor of it for paired-end data. Duplication rate tells you a lot about the original complexity of each sample and potential impact of sequencing depth.

__**However, I would never do PCR duplicate removal on Single-End reads!**__

Many other read de-duplication algorithms rely on mapping position to identify duplicated reads (although some other reference free methods do exist [https://doi.org/10.1186/s12859-016-1192-5](https://doi.org/10.1186/s12859-016-1192-5)). Reads that are mapped to the same position on the genome probably represent the same original fragment sequenced multiple times (think "technical replicates").

However, this approach requires that there be a reference to map reads against and requires that someone maps them!

hts_SuperDeduper does not require a reference or mapped reads. Instead it uses a small portion of each paired read to identify duplicates. If an identical pattern is identified in multiple reads, extra copies are discarded.


<img src="preproc_figures/SD_eval.png" alt="SD_eval" width="80%"/>



<img src="preproc_figures/SD_performance.png" alt="SD_performance" width="80%"/>

We calculated the Youden Index for every combination tested and the point that acquired the highest index value (as compared to Picard MarkDuplicates) occurred at a start position at basepair 5 and a length of 10bp (20bp total over both reads). Though defaults in hts_SuperDeduper are start position at basepair 10 and a length of 10bp.

### Adapter trimming by overlapping reads.

Consider the three scenarios below

**Insert size > length of the number of cycles**

<img src="preproc_figures/overlap_pairs.png" alt="overlap_pairs" width="80%"/>

hts_AdapterTrimmer product: original pairs

hts_Overlapper product: original pairs

**Insert size < length of the number of cycles (10bp min)**

<img src="preproc_figures/overlap_single.png" alt="overlap_single" width="80%"/>

hts_AdapterTrimmer product: original pairs

hts_Overlapper product: extended, single

**Insert size < length of the read length**

<img src="preproc_figures/overlap_adapter.png" alt="overlap_adapter" width="80%"/>

hts_AdapterTrimmer product: adapter trimmed, pairs

hts_Overlapper product: adapter trimmed, single

Both hts_AdapterTrimmer and hts_Overlapper employ this principle to identify and remove adapters for paired-end reads. For paired-end reads the difference between the two are the output, as overlapper produces single-end reads when the pairs overlap and adapter trimmer keeps the paired end format. For single-end reads, adapter trimmer identifies and removes adapters by looking for the adapter sequence, where overlapper just ignores single-end reads (nothing to overlap).


### Now lets see if we can find evidence of Illumina sequencing adapters in our subset.
Remember that Illumina reads must have P5 and P7 adapters and generally look like this (in R1 orientation):

P5---Read1primer---INSERT---IndexReadprimer--index--P7(rc)

This sequence is P7(rc): ATCTCGTATGCCGTCTTCTGCTTG. It should be at the end of any R1 that contains a full-length adapter sequence.

```bash
cd /share/workshop/epigenetics_workshop/$USER/chipseq_example/HTS_testing
zcat JLDY037E.subset_R1.fastq.gz | grep TCTCGTATGCCGTCTTCTGCTTG
```

* *What did you find?*
* *Do you remember how to count the number of instances?*
* *Roughly, what percentage of this data has adapters?*


### Q-window trimming.

As a sequencing run progresses the quality scores tend to get worse. Quality scores are essentially a guess about the accuracy of a base call, so it is common to trim of the worst quality bases.

<img src="preproc_figures/Qwindowtrim.png" alt="Qwindowtrim" width="80%"/>

This is how reads commonly look, they start at "good" quality, increase to "excellent" and degrade to "poor", with R2 always looking worse (except when they don't) than R1 and get worse as the number of cycles increases.

hts_QWindowTrim trims 5' and/or 3' end of the sequence using a windowing (average quality in window) approach.

### What does all this preprocessing get you

Comparing RNAseq mapping count data with raw and preprocessed reads, as an example.

<img src="preproc_figures/final.png" alt="final" width="40%"/>

### Lets put it all together

```bash
cd /share/workshop/epigenetics_workshop/$USER/chipseq_example/HTS_testing

hts_Stats -L JLDY037E_htsStats.json -N "initial stats" \
    -1 JLDY037E.subset_R1.fastq.gz \
    -2 JLDY037E.subset_R2.fastq.gz | \
hts_SeqScreener -A JLDY037E_htsStats.json -N "screen phix" | \
hts_SuperDeduper -A JLDY037E_htsStats.json -N "remove PCR duplicates" | \
hts_AdapterTrimmer -A JLDY037E_htsStats.json -N "trim adapters" | \
hts_NTrimmer -A JLDY037E_htsStats.json -N "remove any remaining 'N' characters" | \
hts_QWindowTrim -A JLDY037E_htsStats.json -N "quality trim the ends of reads" | \
hts_LengthFilter -A JLDY037E_htsStats.json -N "remove reads < 50bp" \
    -n -m 50 | \
hts_Stats -A JLDY037E_htsStats.json -N "final stats" \
    -f JLDY037E.htstream
```

Note the patterns:
* In the first routine we use -1 and -2 to specify the original reads.
* In the final routine -f fastq prefix to write out new preprocessed reads.
* For the log, we specify -L in the first app to write out to a new log, and then use -A for the second routine onward to append log output, generating a single log file at the end.
* All other parameters are algorithm specific, can review using --help

**Questions**
* *Review the final json output, how many reads do we have left?*

* *Confirm that number by counting the number of reads in the final output files.*

* *How many adapters did we detect, cut off?*

* *How many PCR duplicates were there?*

* *Anything else interesting?*

## Run HTStream on the ChIPSeq Project.

We can now run the preprocessing routine across all samples on the real data using a SLURM script, [hts_preproc.slurm](../../software_scripts/scripts/hts_preproc.slurm), that we should take a look at now.

```bash
cd /share/workshop/epigenetics_workshop/$USER/chipseq_example  # We'll run this from the main directory
wget https://ucdavis-bioinformatics-training.github.io/2020-Epigenetics_Workshop/software_scripts/scripts/hts_preproc.slurm -O hts_preproc.slurm
less hts_preproc.slurm
```

When you are done, type "q" to exit.

<div class="script">#!/bin/bash

#SBATCH --job-name=htstream # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=9
#SBATCH --time=12:00:00
#SBATCH --mem=15000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production
#SBATCH --account=epigenetics # cluster account to use for the job
#SBATCH --reservation=epigenetics-workshop# cluster account reservation
#SBATCH --array=1-8
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
</div>

Double check to make sure that slurm_out and 01-HTS_Preproc directories have been created for output, then after looking at the script, let's run it.

```bash
cd /share/workshop/epigenetics_workshop/$USER/chipseq_example
mkdir -p slurm_out  # -p tells mkdir not to complain if the directory already exists
mkdir -p 01-HTS_Preproc
sbatch hts_preproc.slurm  # moment of truth!
```

We can watch the progress of our task array using the 'squeue' command. Takes about 2:30 hours to process each sample.

```bash
squeue -u $USER  # use your username
```

## Now run HTStream on the ATACseq project

Double check to make sure that slurm_out and 01-HTS_Preproc directories have been created for output, then after looking at the script, let's run it.

```bash
cd /share/workshop/epigenetics_workshop/$USER/atacseq_example
wget https://ucdavis-bioinformatics-training.github.io/2020-Epigenetics_Workshop/software_scripts/scripts/hts_preproc.slurm 
-O hts_preproc.slurm

mkdir -p slurm_out  # -p tells mkdir not to complain if the directory already exists
mkdir -p 01-HTS_Preproc
```

**What needs to be changed in the slurm file?**

```bash
sbatch hts_preproc.slurm  # moment of truth!
```

We can watch the progress of our task array using the 'squeue' command. Takes about 1 hour to process each sample.

```bash
squeue -u $USER  # use your username
```

## Quality Assurance - Preprocessing statistics as QA/QC.

Beyond generating "better" data for downstream analysis, cleaning statistics also give you an idea as to the original quality and complexity of the sample, library generation, and sequencing quality.

This can help inform you of how you might change your protocol/procedures in the future, either sample preparation, or in library preparation.

I’ve found it best to perform QA/QC on both the run as a whole (poor samples can affect other samples) and on the samples themselves as they compare to other samples **(BE CONSISTENT!)**.

Reports such as Basespace for Illumina, are great ways to evaluate the run as a whole, the sequencing provider usually does this for you. Plots of the preprocessing summary are a great way to look for technical bias across your experiment. Poor quality samples often appear as outliers and can ethically be removed due to identified technical issues.

1. Let's make sure that all jobs completed successfully.

    Lets first check all the "htstream_%\*.out" and "htstream_%\*.err" files:

    ```bash
    cd /share/workshop/epigenetics_workshop/$USER/chipseq_example
    cat slurm_out/htstream_*.out
    ```

    Look through the output and make sure you don't see any errors. Now do the same for the err files:

    ```bash
    cat slurm_out/htstream_*.err
    ```

    Also, check the output files. First check the number of forward and reverse output files (should be 7 each):

    ```bash
    cd 01-HTS_Preproc
    ls */*R1* | wc -l
    ls */*R2* | wc -l
    ```

    Check the sizes of the files as well. Make sure there are no zero or near-zero size files and also make sure that the size of the files are in the same ballpark as each other:

    ```bash
    ls -lh *
    ```

    **IF** for some reason it didn't finish, is corrupted or you missed the session, please let one of us know and we will help, and you can copy over a completed copy

    ```bash
    #cp -r /share/biocore/workshops/2020_Epigenetics/ChIPseq/HTS_testing /share/workshop/epigenetics_workshop/$USER/chipseq_example/.
    #cp -r /share/biocore/workshops/2020_Epigenetics/ChIPseq/01-HTS_Preproc /share/workshop/epigenetics_workshop/$USER/chipseq_example/.
    ```

1. Let's take a look at the differences in adapter content between the input and output files. First look at the input file:

    ```bash
    cd /share/workshop/epigenetics_workshop/$USER/chipseq_example
    zless 00-RawData/JLDY037E/JLDY037E_S5_L005_R1_001.fastq.gz
    ```

    Let's search for the adapter sequence. Type '/' (a forward slash), and then type **AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC** (the first part of the forward adapter). Press Enter. This will search for the sequence in the file and highlight each time it is found. You can now type "n" to cycle through the places where it is found. When you are done, type "q" to exit. Alternatively, you can use zcat and grep like we did earlier.

    Now look at the output file:

    ```bash
    zless 01-HTS_Preproc/JLDY037E/JLDY037E_R1.fastq.gz
    ```

    If you scroll through the data (using the spacebar), you will see that some of the sequences have been trimmed. Now, try searching for **AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC** again. You shouldn't find it (adapters were trimmed remember), but rarely is anything perfect. You may need to use Control-C to get out of the search and then "q" to exit the 'less' screen.

    Lets grep for the sequence and count occurrences

    ```bash
    zcat  00-RawData/JLDY037E/JLDY037E_S5_L005_R1_001.fastq.gz | grep  AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC | wc -l
    zcat  01-HTS_Preproc/JLDY037E/JLDY037E_R1.fastq.gz | grep  AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC | wc -l
    ```

    * *What is the reduction in adapters found?*

* Perform all the same operations for the ATACseq dataset

If they didn't finish and you need to copy over my copy

```bash
#cp -r /share/biocore/workshops/2020_Epigenetics/ATACseq/01-HTS_Preproc /share/workshop/epigenetics_workshop/$USER/atacseq_example/.
```

1. MultiQC QA/QC Summary of the json files.

Finally lets use [MultiQC](https://multiqc.info/) to generate a summary of our output. Currently MultiQC support for HTStream is in development by Bradley Jenner, and has not been included in the official MultiQC package. If you'd like to try it on your own data, you can find a copy here [https://github.com/s4hts/MultiQC](https://github.com/s4hts/MultiQC).

```bash
## Run multiqc to collect statistics and create a report:
cd /share/workshop/epigenetics_workshop/$USER/chipseq_example
module load multiqc/htstream.dev0
mkdir -p 01-HTS-multiqc-report
multiqc -i ChIPseq-cleaning-report -o 01-HTS-ChIPseq-report ./01-HTS_Preproc
```

**Do the same for the ATACseq experiment**

Transfer ChIPseq-cleaning-report_multiqc_report.html and ATACseq-cleaning-report_multiqc_report.html to your computer and open it in a web browser.

Or in case of emergency, download this copy: [ChIPseq-cleaning-report_multiqc_report.html](ChIPseq-cleaning-report_multiqc_report.html) and [ATACseq-cleaning-report_multiqc_report.html](ATACseq-cleaning-report_multiqc_report.html) for the ATACseq

**Questions**
* *Any problematic samples?*

* *Anything else worth discussing?*

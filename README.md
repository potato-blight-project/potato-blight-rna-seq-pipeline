# potato-blight-rna-seq-pipeline

RNASeq pipeline for the potato blight study.

## Prerequisites

The scripts assume the usage of slurm within the University of Birmingham.

## Reads quality check and trimming

> Software: TrimGalore

This will create a directory called `trimmed` including the trimmed reads:

```
$ sbatch trim_qc.sh
```

To change the running account use:

```
$ sbatch --account=<your account> trim_qc.sh
```

## Align reads

> Software: STAR

The script will align the raw reads against the potato genome using the trimmed reads from the previous step.

The output is a BAM file sorted by coordinates.

```
$ sbatch align_reads.sh
```

## Quantification

> Software: htseq-count

The script will count how many reads map to each feature (intervals on a chromosome).

## DESeq2 analysis

TODO.

# brentlab Novoalign + HTSeq rnaseq pipeline

This runs the 'legacy' brentlabRnaSeq pipeline, which uses the current cluster version of novoalign and htseq. The underlying software is not updated -- samples not part of a current 'experiment set' should not use this pipeline.

## Installation

All you need is an environment with nextflow. This is one way

```{bash}
ml miniconda3

mkdir /scratch/mblab/$USER/conda_envs/nextflow

conda create -p /scratch/mblab/$USER/conda_envs/nextflow nextflow
```

When that finishes, start the environment and update

```{bash}
source activate /scratch/mblab/$USER/conda_envs/nextflow
nextflow self-update
```

## Create a params file

open a text file in your favorite plain text editor.
Paste in the following, with the correct paths.

```{raw}

params {
    output_dir                      = /path/to/output/dir
    sample_sheet                    = /path/to/samplesheet.csv
    KN99_novoalign_index            = /path/to/kn99_novoalign_index
    KN99_stranded_annotation_file   = /path/to/stranded.gff
    KN99_unstranded_annotation_file = /path/to/unstranded.gff
    htseq_count_feature             = exon
}

```

here is an example

```{raw}

params {
    output_dir                      = /scratch/mblab/$USER/rnaseq_pipeline
    sample_sheet                    = /scratch/mblab/$USER/rnaseq_pipeline/sample_sheets/run_12345.csv
    KN99_novoalign_index            = /scratch/mblab/$USER/rnaseq_pipeline/genome_files/KN99/novoalign.index
    KN99_stranded_annotation_file   = /scratch/mblab/$USER/rnaseq_pipeline/genome_files/KN99/stranded.gff
    KN99_unstranded_annotation_file = /scratch/mblab/$USER/rnaseq_pipeline/genome_files/KN99/unstranded.gff
    htseq_count_feature             = exon
}

```

## Run the pipeline

First, copy the data from `lts` to `scratch`. Pay attention to the trailing `/` when using rsync -- they do matter

```{bash}
rsync -aHv `lts/mblab/sequence_data/rnaseq_data/lts_sequence/run_<number>_samples /scratch/mblab/$USER/rnaseq_samples/
```

Next, run the pipeline

```{bash}
nextflow run /path/to/brentlab_rnaseqpipe_nf/main.nf -c params.json
```

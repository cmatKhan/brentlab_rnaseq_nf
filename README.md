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
    output_dir                      = /path/to/output_dir
    sample_sheet                    = /path/to/samplesheet.csv
    KN99_novoalign_index            = /path/to/kn99_novoalign_index
    KN99_fasta                      = /path/to/kn99.fasta
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
    KN99_fasta                      = /scratch/mblab/$USER/rnaseq_pipeline/genome_files/KN99/kn99.fasta
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

To run the pipeline, again open your favorite text editor and create a batch script like so:

```{bash}
#!/bin/bash
#SBATCH --mem=5G
#SBATCH -o rnaseq_pipe_nf.out
#SBATCH -J rna_nf

ml miniconda3 # or whatever conda you use, or however you get to a place where you can use nextflow

source activate /scratch/mblab/$USER/conda_envs/nextflow nextflow

# this is important
mkdir tmp

nextflow run /path/to/brentlab_rnaseqpipe_nf/main.nf -c path/to/your_params.json
```
Submit it

```{bash}
sbatch /path/to/script.sh
```

## Archive the results

You can check progress like so (assuming you're in the same directory from which you launch)

```{bash}
tail -50 rnaseq_pipe_nf.out
```

When the pipeline completes without error, you should move the run_<number>_samples __output__ into `/lts`.

The trailing `/` are important in `rsync`. If you are unsure, ask.

```{bash}
rsync -aHv /path/to/output_dir/rnaseq_pipeline_results/run_<number>_samples lts/mblab/sequence_data/rnaseq_data/lts_align_expr
```
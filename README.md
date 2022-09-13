# brentlab Novoalign + HTSeq rnaseq pipeline

This runs on __only__ the 'new' partition of HTCF
## Installation

All you need is an environment with nextflow. This is one way

```{bash}
# launch interactive session. Note that the memory is just a guess,
# you likely need less. It unlikely that you need more
srun --mem 3GB -n1 --pty bash

# this is a cluster wide installation, or choose whatever you
# typically use
ml miniconda3
```
I recommend that you use your $USER name for your scratch dir 
name as it makes other tasks easier. I also recommend that 
you keep your environments in scratch. [Here is a script/instructions 
to set up a cron job to touch your scratch space every week to ensure 
that the cluster doesn't delete your files](https://github.com/BrentLab/brentlabRnaSeqTools/blob/main/inst/bash/update_timestamps_on_cluster.sh)
```{bash
# assuming 
# 1. your scratch space directory is named your username (see echo $USER)
# 2. you want to put the conda env in a subdir of your scratch called conda_envs
mkdir -p /scratch/mblab/$USER/conda_envs

# create the environment with nextflow in it
conda create -p /scratch/mblab/$USER/conda_envs/nextflow nextflow

# once finished, update nextflow
conda activiate /scratch/mblab/$USER/conda_envs/nextflow
nextflow self-update
```

## Create a params file

open a text file in your favorite plain text editor.
Paste in the following, with the correct paths.

```{json}
{
"output_dir": "/path/to/output_dir",
"sample_sheet": "/path/to/samplesheet.csv",
"run_number": "4689",
"KN99_novoalign_index": "/path/to/kn99_novoalign_index",
"KN99_fasta": "/path/to/kn99.fasta",
"KN99_stranded_annotation_file": "/path/to/stranded.gff",
"KN99_unstranded_annotation_file": "/path/to/unstranded.gff",
"htseq_count_feature": "exon"
}

```

here is an example

```{json}
{
"output_dir": ".",
"sample_sheet": "run_4689.csv",
"run_number": "4689",
"KN99_novoalign_index": "/scratch/mblab/chasem/rnaseq_pipeline/genome_files/KN99/KN99_genome_fungidb.nix",
"KN99_fasta": "/scratch/mblab/chasem/rnaseq_pipeline/genome_files/KN99/KN99_genome_fungidb.fasta",
"KN99_stranded_annotation_file": "/scratch/mblab/chasem/rnaseq_pipeline/genome_files/KN99/KN99_stranded_annotations_fungidb_augment.gff",
"KN99_unstranded_annotation_file": "/scratch/mblab/chasem/rnaseq_pipeline/genome_files/KN99/KN99_no_strand_annotations_fungidb_augment.gff",
"htseq_count_feature": "exon"
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

source activate /scratch/mblab/$USER/conda_envs/nextflow

# this is important, don't change or delete it
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

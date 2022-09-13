# brentlab Novoalign + HTSeq rnaseq pipeline

This runs on __only__ the 'new' partition of HTCF
## Installation

If you are in the brentlab, then you do not need to install any dependencies. 
just pull this respository with `git clone`. If you have already done this, 
then `cd` into the respository and update with `git pull`.
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
    "output_dir": "align_count_results\/run_6019",
    "sample_sheet": "samplesheets\/run_6019.csv",
    "run_number": "6019",
    "KN99_novoalign_index": "\/ref\/mblab\/data\/KN99\/KN99_genome_fungidb.nix",
    "KN99_fasta": "\/ref\/mblab\/data\/KN99\/KN99_genome_fungidb.fasta",
    "KN99_stranded_annotation_file": "\/ref\/mblab\/data\/KN99/KN99_stranded_annotations_fungidb_augment.gff",
    "KN99_unstranded_annotation_file": "\/ref\/mblab\/data\/KN99\/KN99_no_strand_annotations_fungidb_augment.gff",
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

eval $(spack load --sh openjdk)
eval $(spack load --sh nextflow@22.04.5)

# this is important, don't change or delete it
mkdir tmp

nextflow run /path/to/brentlab_rnaseq_nf/main.nf -params-file path/to/your_params.json
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

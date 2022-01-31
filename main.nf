
/* 
 * Legacy Brentlab RnaSeq pipeline, using HTCF software
 * 
 */

/* 
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

log.info """\
Brentlab_legacy_Rnaseq  -  N F    v 2.1 
================================
output_dir                      : $params.output_dir
sample_sheet                    : $params.sample_sheet
KN99_novoalign_index            : $params.KN99_novoalign_index
KN99_fasta                      : $params.KN99_fasta
KN99_stranded_annotation_file   : $params.KN99_stranded_annotation_file
KN99_unstranded_annotation_file : $params.KN99_unstranded_annotation_file
htseq_count_feature             : $params.htseq_count_feature
"""

/* 
 * Import modules 
 */
include { 
  FASTQC;
  NOVOALIGN;
  HTSEQ_EXON;
  HTSEQ_CDS;
  BAM_INDEX;
  MULTIQC } from './modules.nf' 

/* 
 * main pipeline logic
 */
workflow {

    Channel
    .fromPath(params.sample_sheet)
    .splitCsv(header:true)
    .map{row-> tuple(file(row.fastqFilePath), row.strandedness) }
    .set { fastqc_input_ch }
      
      // PART 1: FastQC
      FASTQC(fastqc_input_ch)

      // PART 2: Align
      NOVOALIGN(FASTQC.out[0])
      
      // PART 3a: Count
      HTSEQ_EXON(NOVOALIGN.out[0])

      // PART 3b: Count -- cds
      HTSEQ_CDS(HTSEQ_EXON.OUT[0])

      // PART 4: Index Bams
      BAM_INDEX(HTSEQ_CDS.out[0])
      
      // PART 5: Multiqc
      MULTIQC(FASTQC.out[1].mix(NOVOALIGN.out[1]).mix(HTSEQ_EXON.out[1]).collect())

}

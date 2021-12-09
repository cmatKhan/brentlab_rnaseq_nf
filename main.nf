
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
  HTSEQ;
  BAM_INDEX;
  MULTIQC } from './modules.nf' 

/* 
 * main pipeline logic
 */
workflow {

    Channel
    .fromPath(params.sample_sheet)
    .splitCsv(header:true)
    .map{row-> tuple(row.runNumber, file(row.fastqFilePath), row.strandedness) }
    .set { fastqc_input_ch }
      
      // PART 1: FastQC
      FASTQC(fastqc_input_ch)

      // PART 1: Align
      NOVOALIGN(FASTQC.out[0])
      
      // PART 2: Count
      HTSEQ(NOVOALIGN.out[0])

      // PART 3: Index Bams
      BAM_INDEX(HTSEQ.out[0])
      
      // PART 3: Multiqc
      MULTIQC(FASTQC.OUT[1].mix(NOVOALIGN.OUT[1]).mix(HTSEQ.out[1]).collect())

}

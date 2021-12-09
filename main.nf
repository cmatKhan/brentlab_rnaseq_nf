/*
 * Copyright (c) 2021, Chase Mateusiak.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
 * 
 * This Source Code Form is "Incompatible With Secondary Licenses", as
 * defined by the Mozilla Public License, v. 2.0.
 */


/* 
 * bartNPNF is a proof concept pipeline for nextflow
 * 
 * Chase Mateusiak
 */

/* 
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

log.info """\
npNF  -  N F    v 2.1 
================================
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
    .fromPath(params.fastq_file_list)
    .splitCsv(header:true)
    .map{row-> tuple(row.runNumber, file(row.fastqFilePath), row.strandedness) }
    .set { fastqc_input_ch }

      
      // PART 1: Align
      NOVOALIGN(genes_ch,
                    params.nps,
                    params.ntree, 
                    params.test_data)
      
      // PART 2: Count
      HTSEQ(genes_ch,
                    params.nps,
                    params.ntree, 
                    params.test_data)
      
      // PART 3: Multiqc
      MULTIQC(genes_ch,
                    params.nps,
                    params.ntree, 
                    params.test_data)

}

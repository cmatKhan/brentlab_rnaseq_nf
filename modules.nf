/*
 * Process 1: FASTQC
 */

process FASTQC {
   
    beforeScript "ml fastqc/0.11.7-java-11"

    input:
    tuple val(runNumber), path(fastqFilePath), val(strandedness) from fastqc_input_ch

    output:
    tuple val(runNumber), path(fastqFilePath), val(strandedness) from fastq_align_ch
    file("${fastqc_out}") into fastqc_ch


    script:
    fastqc_out="${params.output_dir}/rnaseq_pipeline_results/run_${runNumber}_samples/{align,logs,count,fastqc}"
    """
    mkdir -p ${fastqc_out}
    fastqc -o ${fastqc_out} -f fastq -q ${fastqFilePath}
    """  
}  

/*
 * Process 1: FASTQC
 */

process NOVOALIGN {

    label 'align_count'

    beforeScript "ml novoalign/3.09.01 samtools"
    publishDir "${params.output_dir}/rnaseq_pipeline_results/run_${runNumber}_samples/logs", overwite: true, pattern: "*.log"


    input:
        tuple val(runNumber), path(fastqFilePath), val(strandedness) from fastq_align_ch
    output:
        tuple val(runNumber), val(strandedness), val(fastq_simple_name), path("*_sorted_aligned_reads.bam") into bam_align_ch
        path "*_novosort.log" novoalign_log_ch

    script:
        fastq_simple_name = fastqFilePath.getSimpleName()
            """
            novoalign -r All \\
                      -c 8 \\
                      -o SAM \\
                      -d ${params.KN99_novoalign_index} \\
                      -f ${fastqFilePath} 2> ${fastq_simple_name}_novoalign.log | \\
            samtools view -bS | \\
            novosort - \\
                     --threads 8 \\
                     --markDuplicates \\
                     --index \\
                     -o ${fastq_simple_name}_sorted_aligned_reads.bam 2> ${fastq_simple_name}_novosort.log
            """
}

/*
 * Process 2: HTSEQ
 */

process HTSEQ {

    
    label 'align_count'

    beforeScript "ml samtools htseq/0.9.1"

    publishDir "${params.output_dir}/rnaseq_pipeline_results/run_${runNumber}_samples/logs", overwite: true, pattern: "*.log"
    publishDir "${params.output_dir}/rnaseq_pipeline_results/run_${runNumber}_samples/count", overwite: true, pattern: "*.tsv"
    publishDir "${params.output_dir}/rnaseq_pipeline_results/run_${runNumber}_samples/align", overwite: true, pattern: "*_sorted_aligned_reads_with_annote.bam"


    input:
        tuple val(runNumber), val(strandedness), val(fastq_simple_name), path(sorted_bam) from bam_align_ch
    output:
        path "*_sorted_aligned_reads_with_annote.bam" into bam_align_with_htseq_annote_ch
        path "*.tsv" into htseq_ch
        path "*.log"

    script:

        if (strandedness == 'reverse')
            """
            htseq-count -f bam \\
                        -o ${fastq_simple_name}_htseq_annote.sam \\
                        -s ${strandedness} \\
                        -t ${params.htseq_count_feature} \\
                        -i gene \\
                        ${sorted_bam} \\
                        ${params.KN99_stranded_annotation_file} \\
                        1> ${fastq_simple_name}_read_count.tsv 2> ${fastq_simple_name}_htseq.log

            sed "s/\t//" ${fastq_simple_name}_htseq_annote.sam > ${fastq_simple_name}_no_tab_sam.sam

            samtools view --threads 8 ${sorted_bam} | \\
            paste - ${fastq_simple_name}_no_tab_sam.sam | \\
            samtools view --threads 8 -bS -T ${params.KN99_genome} > ${fastq_simple_name}_sorted_aligned_reads_with_annote.bam

            """
        else if (strandedness == 'unstranded')
            """
            htseq-count -f bam \\
                        -o ${fastq_simple_name}_htseq_annote.sam \\
                        -s ${strandedness} \\
                        -t ${params.htseq_count_feature} \\
                        -i gene \\
                        ${sorted_bam} \\
                        ${params.KN99_unstranded_annotation_file} \\
                        1> ${fastq_simple_name}_read_count.tsv 2> ${fastq_simple_name}_htseq.log

            sed "s/\t//" ${fastq_simple_name}_htseq_annote.sam > ${fastq_simple_name}_no_tab_sam.sam

            samtools view --threads 8 ${sorted_bam} | \\
            paste - ${fastq_simple_name}_no_tab_sam.sam | \\
            samtools view --threads 8 -bS -T ${params.KN99_genome} > ${fastq_simple_name}_sorted_aligned_reads_with_annote.bam

            """
}

/*
 * Process 3: INDEX
 */


process BAM_INDEX {


    label 'index'
    beforeScript "ml novoalign/3.09.01 samtools"
    publishDir "${params.output_dir}/rnaseq_pipeline_results/run_${runNumber}_samples/align", overwite: true, pattern: "*.bai"


    input:
        path(bam) from bam_align_with_htseq_annote_ch
    output:
        path "*.bai"

    script:

            """
            samtools index $bam
            """
}

/*
 * Process 4: MULTIQC
 */

process MULTIQC {

    label 'multiqc'

    publishDir "${params.output_dir}/rnaseq_pipeline_results/run_${runNumber}_samples", pattern:"*.html"
       
    input:
    file('*') from novoalign_log_ch.mix(htseq_ch).mix(fastqc_ch).collect()
    
    output:
    file('multiqc_report.html')  
     
    script:
    """
    multiqc . 
    """
}


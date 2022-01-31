/*
 * Process 1: FASTQC
 */

process FASTQC {

    label 'fastqc'
   
    beforeScript "ml fastqc/0.11.7-java-11"

    input:
    tuple path(fastqFilePath), val(strandedness)

    output:
    tuple path(fastqFilePath), val(strandedness)
    file("${fastqc_out}/*") 


    script:
    fastqc_out="fastqc_out"
    """
    mkdir -p ${fastqc_out}
    fastqc -o ${fastqc_out} -t 8 -f fastq -q ${fastqFilePath}
    """  
}  

/*
 * Process 2: NOVOALIGN
 */

process NOVOALIGN {

    label 'align'

    beforeScript "ml novoalign/3.09.01 samtools"
    publishDir "${params.output_dir}/rnaseq_pipeline_results/run_${params.run_number}_samples/logs", overwite: true, pattern: "*.log", mode: 'copy'


    input:
        tuple path(fastqFilePath), val(strandedness)
    output:
        tuple val(strandedness), val(fastq_simple_name), path("*_sorted_aligned_reads.bam") 
        path "*.log" 

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
 * Process 3a: HTSEQ_EXON
 */

process HTSEQ_EXON {

    
    label 'htseq'

    beforeScript "ml samtools htseq/0.9.1"

    publishDir "${params.output_dir}/rnaseq_pipeline_results/run_${params.run_number}_samples/logs", overwite: true, pattern: "*.log", mode: 'copy'
    publishDir "${params.output_dir}/rnaseq_pipeline_results/run_${params.run_number}_samples/count", overwite: true, pattern: "*.tsv", mode: 'copy'


    input:
        tuple val(strandedness), val(fastq_simple_name), path(sorted_bam)
    output:
        tuple val(strandedness), val(fastq_simple_name), path("*_sorted_aligned_reads_with_annote.bam")
        path "*.tsv" 
        path "*.log"

    script:

        if (strandedness == 'reverse')
            """
            htseq-count -f bam \\
                        -o ${fastq_simple_name}_htseq_annote.sam \\
                        -s ${strandedness} \\
                        -t exon \\
                        -i gene \\
                        ${sorted_bam} \\
                        ${params.KN99_stranded_annotation_file} \\
                        1> ${fastq_simple_name}_read_count.tsv 2> ${fastq_simple_name}_htseq.log

            sed "s/\t//" ${fastq_simple_name}_htseq_annote.sam > ${fastq_simple_name}_no_tab_sam.sam

            samtools view --threads 8 ${sorted_bam} | \\
            paste - ${fastq_simple_name}_no_tab_sam.sam | \\
            samtools view --threads 8 -bS -T ${params.KN99_fasta} > \\
            ${fastq_simple_name}_sorted_aligned_reads_with_annote.bam

            """
        else 
            """
            htseq-count -f bam \\
                        -o ${fastq_simple_name}_htseq_annote.sam \\
                        -s no \\
                        -t ${params.htseq_count_feature} \\
                        -i gene \\
                        ${sorted_bam} \\
                        ${params.KN99_unstranded_annotation_file} \\
                        1> ${fastq_simple_name}_read_count.tsv 2> ${fastq_simple_name}_htseq.log

            sed "s/\t//" ${fastq_simple_name}_htseq_annote.sam > ${fastq_simple_name}_no_tab_sam.sam

            samtools view --threads 8 ${sorted_bam} | \\
            paste - ${fastq_simple_name}_no_tab_sam.sam | \\
            samtools view --threads 8 -bS -T ${params.KN99_fasta} > \\
            ${fastq_simple_name}_sorted_aligned_reads_with_annote.bam

            """
}


/*
 * Process 3b: HTSEQ_CDS
 */

process HTSEQ_CDS {

    
    label 'htseq'

    beforeScript "ml samtools htseq/0.9.1"

    publishDir "${params.output_dir}/rnaseq_pipeline_results/run_${params.run_number}_samples/logs", overwite: true, pattern: "*.log", mode: 'copy'
    publishDir "${params.output_dir}/rnaseq_pipeline_results/run_${params.run_number}_samples/count", overwite: true, pattern: "*.tsv", mode: 'copy'


    input:
        tuple val(strandedness), val(fastq_simple_name), path(sorted_bam)
    output:
        path(sorted_bam)
        path "*.tsv" 
        path "*.log"

    script:

        if (strandedness == 'reverse')
            """
            htseq-count -f bam \\
                        -o ${fastq_simple_name}_htseq_annote.sam \\
                        -s ${strandedness} \\
                        -t cds \\
                        -i gene \\
                        ${sorted_bam} \\
                        ${params.KN99_stranded_annotation_file} \\
                        1> ${fastq_simple_name}_read_count_cds.tsv 2> ${fastq_simple_name}_htseq_cds.log
            """
        else 
            """
            htseq-count -f bam \\
                        -o ${fastq_simple_name}_htseq_annote.sam \\
                        -s no \\
                        -t ${params.htseq_count_feature} \\
                        -i gene \\
                        ${sorted_bam} \\
                        ${params.KN99_unstranded_annotation_file} \\
                        1> ${fastq_simple_name}_read_count.tsv 2> ${fastq_simple_name}_htseq.log
            """
}

/*
 * Process 4: INDEX
 */


process BAM_INDEX {


    label 'index'
    beforeScript "ml novoalign/3.09.01 samtools"

    publishDir "${params.output_dir}/rnaseq_pipeline_results/run_${params.run_number}_samples/align", overwite: true, pattern: "*.bam*", mode: 'copy'

    input:
        tuple path(bam)
    output:
        path(bam)
        path "*.bai*"

    script:

            """
            samtools index -@ 8 $bam
            """
}

/*
 * Process 5: MULTIQC
 */

process MULTIQC {

    label 'multiqc'

    beforeScript "ml multiqc"

    publishDir "${params.output_dir}/rnaseq_pipeline_results/run_${params.run_number}_samples/multiqc", pattern:"*", mode: 'copy'
       
    input:
    file('*')
    
    output:
    path "*multiqc_report.html"
    path "*_data"       
     
    script:
    """
    multiqc . 
    """
}



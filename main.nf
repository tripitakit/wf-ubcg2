#!/usr/bin/env nextflow

// wf-ubcg2
//
// Combined workflow for bacterial phylogenomics:
// - Stage 1: Extract Universal Bacterial Core Genes (UBCG2) from genome assemblies
// - Stage 2: Align and concatenate UCG genes into phylogenetic supermatrix
//
// Both stages are independent and can be run separately or together.

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fasta_ingress } from './lib/ingress'
include { getParams } from './lib/common'

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

// ============================================================================
// STAGE 1: UCG Extraction Processes
// ============================================================================

process getVersions_extraction {
    label "wftemplate"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "versions_extraction.txt"
    cpus 1
    output:
        path "versions_extraction.txt"
    script:
    """
    echo "wf-ubcg2,v1.0.0" > versions_extraction.txt
    echo "UBCG2,2.0" >> versions_extraction.txt
    java -version 2>&1 | head -n 1 | sed 's/^/Java,/' >> versions_extraction.txt
    """
}

process extractUCG {
    label "wfucg"
    publishDir "${params.out_dir}/ucg_output", mode: 'copy', pattern: "${meta.alias}/*.ucg"
    cpus 4
    memory '8GB'
    input:
        tuple val(meta), path(fasta)
    output:
        tuple val(meta), path("${meta.alias}/*.ucg"), emit: ucg_file
        path "${meta.alias}/ubcg2.log", optional: true, emit: log
    script:
    """
    # Unset JAVA_TOOL_OPTIONS for Java 8 compatibility
    unset JAVA_TOOL_OPTIONS

    mkdir -p ${meta.alias}

    java -jar /programs/UBCG2.jar \\
        -i ${fasta} \\
        -ucg_dir ${meta.alias} \\
        -label ${meta.alias} \\
        -hmm /programs/hmm/ubcg_v2.hmm \\
        -t ${task.cpus} \\
        2>&1 | tee ${meta.alias}/ubcg2.log
    """
}

// ============================================================================
// STAGE 2: Alignment and Concatenation Processes
// ============================================================================

process getVersions_alignment {
    label "wfucgalign"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "versions_alignment.txt"
    cpus 1
    output:
        path "versions_alignment.txt"
    script:
    """
    echo "wf-ubcg2,v1.0.0" > versions_alignment.txt
    mafft --version 2>&1 | head -n 1 | sed 's/^/MAFFT,/' >> versions_alignment.txt
    python3 --version 2>&1 | sed 's/Python /Python,/' >> versions_alignment.txt
    """
}

process parseUCG {
    label "wfucgalign"
    publishDir path: "${params.out_dir}/parsed_genes", mode: 'copy', pattern: "*.fasta"
    publishDir path: "${params.out_dir}/parsed_genes", mode: 'copy', pattern: "gene_selection_report.txt"
    cpus 1
    memory '2GB'
    input:
        path ucg_files
    output:
        path "*.fasta", emit: gene_fastas
        path "gene_selection_report.txt", emit: report
    script:
    """
    parse_ucg.py . ${ucg_files}
    """
}

process alignGene {
    label "wfucgalign"
    publishDir "${params.out_dir}/alignments", mode: 'copy', pattern: "*.aln"
    cpus params.mafft_threads
    memory '4GB'
    input:
        path gene_fasta
    output:
        path "*.aln", emit: alignment
    script:
    def gene_name = gene_fasta.baseName
    """
    mafft --thread ${task.cpus} --auto ${gene_fasta} > ${gene_name}.aln
    """
}

process trimAlignment {
    label "wfucgalign"
    publishDir "${params.out_dir}/trimmed_alignments", mode: 'copy', pattern: "*.trimmed.fasta"
    cpus 1
    memory '2GB'
    input:
        path alignment
    output:
        path "*.trimmed.fasta", emit: trimmed
    script:
    def gene_name = alignment.baseName
    """
    trim_alignment.py ${alignment} ${gene_name}.trimmed.fasta
    """
}

process concatenateAlignments {
    label "wfucgalign"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "concatenated_supermatrix.*"
    cpus 1
    memory '4GB'
    input:
        path trimmed_alignments
    output:
        path "concatenated_supermatrix.fasta", emit: supermatrix
        path "concatenated_supermatrix.partitions.txt", emit: partitions
    script:
    """
    concatenate_alignments.py concatenated_supermatrix ${trimmed_alignments}
    """
}

// ============================================================================
// Workflow Modules
// ============================================================================

workflow stage1_extraction {
    take:
        fasta_input  // [meta, fasta_path]
    main:
        extractUCG(fasta_input)
        getVersions_extraction()
    emit:
        ucg_file = extractUCG.out.ucg_file
        ucg_log = extractUCG.out.log
        versions = getVersions_extraction.out
}

workflow stage2_alignment {
    take:
        ucg_files_ch  // channel of UCG files
    main:
        parseUCG(ucg_files_ch)
        alignGene(parseUCG.out.gene_fastas.flatten())
        trimAlignment(alignGene.out.alignment)
        concatenateAlignments(trimAlignment.out.trimmed.collect())
        getVersions_alignment()
    emit:
        supermatrix = concatenateAlignments.out.supermatrix
        partitions = concatenateAlignments.out.partitions
        versions = getVersions_alignment.out
}

// ============================================================================
// Main Entrypoint Workflow
// ============================================================================

WorkflowMain.initialise(workflow, params, log)
workflow {

    Pinguscript.ping_start(nextflow, workflow, params)

    // Check that at least one stage is requested
    def run_stage1 = params.fasta != null
    def run_stage2 = params.ucg_dir != null

    if (!run_stage1 && !run_stage2) {
        error """
        ==========================================
        ERROR: No workflow stage specified!

        Please provide at least one of:
          --fasta       Path to FASTA file (for UCG extraction)
          --ucg_dir     Path to directory with UCG files (for alignment)

        Examples:
          # Run only Stage 1 (UCG extraction)
          nextflow run main.nf --fasta genome.fasta

          # Run only Stage 2 (alignment)
          nextflow run main.nf --ucg_dir ucg_output/

          # Run both stages sequentially
          nextflow run main.nf --fasta genome.fasta --ucg_dir ucg_output/
        ==========================================
        """
    }

    // ========================================
    // Stage 1: UCG Extraction
    // ========================================
    if (run_stage1) {
        log.info """
        ==========================================
        Stage 1: UCG Extraction
        ==========================================
        Input FASTA: ${params.fasta}
        Output directory: ${params.out_dir}/ucg_output
        ==========================================
        """

        fasta_data = fasta_ingress([
            "input": params.fasta
        ])

        stage1_extraction(fasta_data)
    }

    // ========================================
    // Stage 2: Alignment and Concatenation
    // ========================================
    if (run_stage2) {
        log.info """
        ==========================================
        Stage 2: Alignment & Concatenation
        ==========================================
        UCG directory: ${params.ucg_dir}
        MAFFT threads: ${params.mafft_threads}
        Output directory: ${params.out_dir}
        ==========================================
        """

        // Find all UCG files in the directory
        Channel
            .fromPath("${params.ucg_dir}/*.ucg", checkIfExists: false)
            .collect()
            .map { files ->
                if (files.size() == 0) {
                    error "No UCG files (*.ucg) found in ${params.ucg_dir}"
                }
                log.info "Found ${files.size()} UCG files:"
                files.each { log.info "  - ${it.name}" }
                return files
            }
            .set { ucg_files_ch }

        stage2_alignment(ucg_files_ch)
    }

    // Get workflow parameters
    getParams()
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)

    log.info """
    ==========================================
    Pipeline Completed!
    ==========================================
    Status: ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Output directory: ${params.out_dir}
    Duration: ${workflow.duration}
    ==========================================
    """
}

workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
    log.error "Pipeline execution failed: ${workflow.errorMessage}"
}

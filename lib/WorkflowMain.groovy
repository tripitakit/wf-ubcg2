// Simplified WorkflowMain for wf-ubcg2

class WorkflowMain {

    // Citation string for pipeline
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
            "* The nf-core framework\n" +
            "  https://doi.org/10.1038/s41587-020-0439-x\n\n"
    }

    // Generate help string
    public static String help(workflow, params, log) {
        String help_string = ''
        help_string += NfcoreTemplate.logo(workflow, params.monochrome_logs)
        help_string += "\nTypical pipeline command:\n\n"
        help_string += "  nextflow run ${workflow.manifest.name} \\\n"
        params.wf.example_cmd.each { cmd ->
            if (cmd && cmd.trim()) {
                help_string += "    ${cmd} \\\n"
            }
        }
        help_string += '\n' + citation(workflow) + '\n'
        return help_string
    }

    // Generate parameter summary log string
    public static String paramsSummaryLog(workflow, params, log) {
        String workflow_version = NfcoreTemplate.version(workflow)
        String summary_log = ''
        summary_log += NfcoreTemplate.logo(workflow, params.monochrome_logs)
        summary_log += NfcoreTemplate.dashedLine(params.monochrome_logs)
        summary_log += "\nWorkflow: ${workflow.manifest.name} ${workflow_version}\n"
        summary_log += "Run Name: ${workflow.runName}\n"
        summary_log += NfcoreTemplate.dashedLine(params.monochrome_logs)
        summary_log += "\nParameters:\n"

        // Show main parameters
        if (params.fasta) summary_log += "  --fasta: ${params.fasta}\n"
        if (params.ucg_dir) summary_log += "  --ucg_dir: ${params.ucg_dir}\n"
        summary_log += "  --out_dir: ${params.out_dir}\n"
        if (params.mafft_threads) summary_log += "  --mafft_threads: ${params.mafft_threads}\n"

        summary_log += NfcoreTemplate.dashedLine(params.monochrome_logs)
        summary_log += '\n' + citation(workflow) + '\n'
        summary_log += NfcoreTemplate.dashedLine(params.monochrome_logs)
        return summary_log
    }

    // Validate parameters and print summary to screen
    public static void initialise(workflow, params, log) {
        // Print help to screen if required
        if (params.help) {
            log.info help(workflow, params, log)
            System.exit(0)
        }

        // Print workflow version and exit on --version
        if (params.version) {
            String workflow_version = NfcoreTemplate.version(workflow)
            log.info "${workflow.manifest.name} ${workflow_version}"
            System.exit(0)
        }

        // Explode on conda
        try {
            if (workflow.session.config.conda.enabled) {
                log.error "Sorry, this workflow is not compatible with Conda, please use -profile standard (Docker) or -profile singularity."
                System.exit(1)
            }
        } catch(Exception e) {}

        // Print parameter summary log to screen
        log.info paramsSummaryLog(workflow, params, log)
    }
}

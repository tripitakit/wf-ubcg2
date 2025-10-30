// Simplified NfcoreTemplate for wf-ubcg2

class NfcoreTemplate {

    // Return workflow version string
    public static String version(workflow) {
        return workflow.manifest.version ?: "unknown"
    }

    // Generate workflow logo
    public static String logo(workflow, monochrome_logs) {
        Map colors = logColours(monochrome_logs)
        String logo = """\n
        ${colors.cyan}=====================================${colors.reset}
        ${colors.blue}  ${workflow.manifest.name}${colors.reset}
        ${colors.cyan}=====================================${colors.reset}
        """.stripIndent()
        return logo
    }

    // Generate a dashed line
    public static String dashedLine(monochrome_logs) {
        Map colors = logColours(monochrome_logs)
        return "${colors.dim}----------------------------------------------------${colors.reset}\n"
    }

    // Define ANSI color codes
    private static Map logColours(Boolean monochrome_logs) {
        Map colorcodes = [:]

        if (monochrome_logs) {
            colorcodes = [
                'reset': '',
                'bold': '',
                'dim': '',
                'underlined': '',
                'blink': '',
                'reverse': '',
                'hidden': '',
                'black': '',
                'red': '',
                'green': '',
                'yellow': '',
                'blue': '',
                'purple': '',
                'cyan': '',
                'white': ''
            ]
        } else {
            colorcodes = [
                'reset': '\033[0m',
                'bold': '\033[1m',
                'dim': '\033[2m',
                'underlined': '\033[4m',
                'blink': '\033[5m',
                'reverse': '\033[7m',
                'hidden': '\033[8m',
                'black': '\033[0;30m',
                'red': '\033[0;31m',
                'green': '\033[0;32m',
                'yellow': '\033[0;33m',
                'blue': '\033[0;34m',
                'purple': '\033[0;35m',
                'cyan': '\033[0;36m',
                'white': '\033[0;37m'
            ]
        }
        return colorcodes
    }
}

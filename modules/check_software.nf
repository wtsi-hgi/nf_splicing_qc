def check_software_exists(tool) {
    try {
        def process = ["which", tool].execute()
        process.waitFor()
        return process.exitValue() == 0
    } catch (Exception e) {
        return false
    }
}

def check_required(required_tools) {
    def missing_tools = required_tools.findAll { !check_software_exists(it) }

    println "*-----------------------------*"
    println "Checking software:"
    required_tools.each { tool ->
        if (check_software_exists(tool)) {
            println "    |----> ${tool} is available"
        } else {
            println "    |----> ${tool} is not found"
        }
    }

    if (missing_tools) {
        println "Error: the following tools are missing: ${missing_tools.join(', ')}"
        System.exit(1)
    }

    println "Done: all required tools are available. Proceeding with the pipeline."
}

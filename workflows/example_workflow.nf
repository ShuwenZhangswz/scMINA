#!/usr/bin/env nextflow

/*
 * Example Nextflow workflow for scMINA
 * Demonstrates Python and R integration
 */

// Enable DSL2
nextflow.enable.dsl=2

// Workflow parameters (DSL2 style)
params.input_data = "/group/gquongrp/workspaces/hongruhu/sea-ad/software/scMINA/data/example_data.h5ad"
params.output_dir = "./results"
params.python_script = "/group/gquongrp/workspaces/hongruhu/sea-ad/software/scMINA/scripts/example_process_data.py"
params.r_script = "/group/gquongrp/workspaces/hongruhu/sea-ad/software/scMINA/scripts/example_visualize_data.R"

// Process 1: Python data processing
process python_processing {
    label 'high_memory'
    
    input:
    path input_file
    path python_script
    
    output:
    path "processed_data.h5ad", emit: processed_data
    path "python_log.txt", emit: log
    
    script:
    """
    python ${python_script} \
        --input ${input_file} \
        --output processed_data.h5ad \
        --log python_log.txt
    """
}

// Process 2: R visualization
process r_visualization {
    label 'high_memory'
    
    input:
    path processed_data
    path r_script
    
    output:
    path "*.png", emit: plots
    path "r_log.txt", emit: log
    
    script:
    """
    Rscript ${r_script} \
        --input ${processed_data} \
        --output . \
        --log r_log.txt
    """
}

// Process 3: Generate report (Python)
process generate_report {
    label 'high_memory'
    
    input:
    path processed_data
    path plots
    path python_log
    path r_log
    
    output:
    path "report.html", emit: report
    
    script:
    '''
    PLOTS_COUNT=$(ls -1 *.png 2>/dev/null | wc -l)
    cat > report.html << EOF
    <html>
    <head><title>scMINA Analysis Report</title></head>
    <body>
    <h1>scMINA Analysis Report</h1>
    <h2>Input Data</h2>
    <p>Processed data: !{processed_data}</p>
    <h2>Results</h2>
    <p>Generated plots: ${PLOTS_COUNT} files</p>
    <h2>Logs</h2>
    <p>Python log: <a href="!{python_log}">$(basename !{python_log})</a></p>
    <p>R log: <a href="!{r_log}">$(basename !{r_log})</a></p>
    </body>
    </html>
    EOF
    '''
}

// Main workflow
workflow {
    // Check if input file exists
    if (!file(params.input_data).exists()) {
        exit 1, "Input file not found: ${params.input_data}"
    }
    
    // Run Python processing
    python_processing(params.input_data, params.python_script)
    
    // Run R visualization
    r_visualization(python_processing.out.processed_data, params.r_script)
    
    // Generate final report
    generate_report(
        python_processing.out.processed_data,
        r_visualization.out.plots,
        python_processing.out.log,
        r_visualization.out.log
    )
    
    // Print completion message
    println """
    ==========================================
    Workflow completed successfully!
    ==========================================
    Results saved to: ${params.output_dir}
    Report: ${generate_report.out.report}
    ==========================================
    """
}

// Handle workflow completion
workflow.onComplete {
    println "Workflow completed at: ${workflow.complete}"
    println "Execution time: ${workflow.duration}"
}

// Handle workflow error
workflow.onError {
    println "Workflow failed: ${workflow.errorMessage}"
    System.exit(1)
}

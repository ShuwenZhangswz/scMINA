#!/usr/bin/env nextflow

/*
 * Example Nextflow workflow for scMINA
 * Demonstrates Python and R integration
 */

// Workflow parameters
params {
    input_data = "data/example.h5ad"
    output_dir = "results"
    python_script = "scripts/example_process_data.py"
    r_script = "scripts/example_visualize_data.R"
}

// Print workflow information
workflow {
    println """
    ==========================================
    scMINA Nextflow Workflow
    ==========================================
    Input data: ${params.input_data}
    Output directory: ${params.output_dir}
    Python script: ${params.python_script}
    R script: ${params.r_script}
    ==========================================
    """
}

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
    """
    python -c "
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Create simple HTML report
html_content = f'''
<html>
<head><title>scMINA Analysis Report</title></head>
<body>
<h1>scMINA Analysis Report</h1>
<h2>Input Data</h2>
<p>Processed data: {processed_data}</p>
<h2>Results</h2>
<p>Generated plots: {len(list(Path('.').glob('*.png')))} files</p>
<h2>Logs</h2>
<pre>Python log:\n{open(python_log).read()}</pre>
<pre>R log:\n{open(r_log).read()}</pre>
</body>
</html>
'''
with open('report.html', 'w') as f:
    f.write(html_content)
"
    """
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

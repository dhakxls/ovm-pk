"""HTML report generation for pipeline results."""
from pathlib import Path
import json

def generate_html_report(results_path: Path, output_path: Path):
    """Create HTML report from pipeline results."""
    with open(results_path) as f:
        results = json.load(f)
    
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>OVM-PK Run Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; line-height: 1.6; }}
        .container {{ max-width: 800px; margin: 0 auto; }}
        .section {{ margin-bottom: 2rem; }}
        .success {{ color: green; }}
        .failed {{ color: red; }}
        pre {{ background: #f4f4f4; padding: 1rem; overflow-x: auto; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>OVM-PK Run Report</h1>
        
        <div class="section">
            <h2>Status</h2>
            <p class="{'success' if results['status'] == 'success' else 'failed'}">
                {results['status'].upper()}
            </p>
            <p>Completed at: {results['timestamp']}</p>
            <p>Execution time: {results.get('performance', {}).get('execution_time_sec', 'N/A')} seconds</p>
        </div>
        
        <div class="section">
            <h2>System Info</h2>
            <pre>{json.dumps(results.get('system_info', {}), indent=2)}</pre>
        </div>
        
        <div class="section">
            <h2>Configuration</h2>
            <pre>{json.dumps(results.get('config', {}), indent=2)}</pre>
        </div>
    </div>
</body>
</html>
"""
    
    output_path.write_text(html_content)

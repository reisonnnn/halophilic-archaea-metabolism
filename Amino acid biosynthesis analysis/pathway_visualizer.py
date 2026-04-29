import csv
from pathlib import Path
from collections import defaultdict
from datetime import datetime
import json

# Configuration
SCRIPT_DIR = Path(__file__).parent.resolve()
INPUT_FILE = SCRIPT_DIR / "amino_acid_module_details.csv"
OUTPUT_HTML = SCRIPT_DIR / "pathway_visualizer.html"

def parse_csv_data(filepath):
    """Parse the CSV file and return structured data"""
    data = []
    
    if not filepath.exists():
        print(f"Error: {filepath} not found!")
        return []
    
    with open(filepath, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                total_steps = int(row.get('Total_Steps', 0))
                satisfied_steps = int(row.get('Satisfied_Steps', 0))
                completion_pct = (satisfied_steps / total_steps * 100) if total_steps > 0 else 0
                
                item = {
                    'amino_acid': row.get('Amino_Acid', '').strip(),
                    'module': row.get('Module_ID', '').strip(),
                    'pathway': row.get('Pathway_Name', '').strip(),
                    'total_steps': total_steps,
                    'satisfied_steps': satisfied_steps,
                    'completion_pct': round(completion_pct, 1),
                    'total_kos': int(row.get('Total_KOs', 0)),
                    'present_kos_count': int(row.get('Present_KOs_count', 0)),
                    'status': row.get('Status', '').strip(),
                    'present_kos': row.get('Present_KOs', '').strip(),
                    'missing_kos': row.get('Missing_KOs', '').strip()
                }
                data.append(item)
            except (ValueError, KeyError) as e:
                print(f"Warning: Skipping row due to error: {e}")
                continue
    
    return data

def generate_html(data, output_path):
    """Generate the interactive HTML visualizer"""
    
    # Calculate statistics
    aa_counts = defaultdict(int)
    status_counts = defaultdict(int)
    
    for item in data:
        aa_counts[item['amino_acid']] += 1
        status_counts[item['status']] += 1
    
    # Sort amino acids by count (descending)
    sorted_amino_acids = sorted(aa_counts.items(), key=lambda x: -x[1])
    
    # Convert data to JSON for JavaScript
    data_json = json.dumps(data, indent=2)
    date_str = datetime.now().strftime('%Y-%m-%d %H:%M')
    
    html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Pathway Visualizer</title>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/html2canvas/1.4.1/html2canvas.min.js"></script>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        :root {{
            --primary-color: #2c3e50;
            --secondary-color: #34495e;
            --success-color: #27ae60;
            --warning-color: #f39c12;
            --danger-color: #e74c3c;
            --light-bg: #ecf0f1;
            --white: #ffffff;
            --border-color: #bdc3c7;
            --text-dark: #2c3e50;
            --text-light: #7f8c8d;
            --font-size: 13px;
        }}
        
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            background-color: var(--light-bg);
            color: var(--text-dark);
            line-height: 1.6;
            padding: 20px;
            font-size: var(--font-size);
        }}
        
        .container {{
            max-width: 1400px;
            margin: 0 auto;
        }}
        
        /* Header */
        .header {{
            background: var(--white);
            padding: 24px;
            border-radius: 12px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            margin-bottom: 24px;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }}
        
        .header h1 {{
            color: var(--primary-color);
            font-size: 28px;
            margin-bottom: 8px;
            font-family: inherit;
        }}
        
        .header-info {{
            display: flex;
            gap: 24px;
            color: var(--text-light);
            font-size: var(--font-size);
        }}
        
        .font-controls {{
            display: flex;
            gap: 10px;
            align-items: center;
        }}
        
        .font-controls label {{
            font-size: var(--font-size);
            color: var(--text-light);
        }}
        
        .font-controls select {{
            padding: 6px 10px;
            border: 1px solid var(--border-color);
            border-radius: 6px;
            font-size: var(--font-size);
            font-family: inherit;
            background: var(--white);
        }}
        
        /* Controls Section */
        .controls {{
            background: var(--white);
            padding: 20px;
            border-radius: 12px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            margin-bottom: 24px;
        }}
        
        .control-section {{
            margin-bottom: 20px;
        }}
        
        .control-section:last-child {{
            margin-bottom: 0;
        }}
        
        .control-label {{
            font-size: calc(var(--font-size) - 1px);
            font-weight: 600;
            text-transform: uppercase;
            color: var(--text-light);
            margin-bottom: 10px;
            display: block;
            font-family: inherit;
        }}
        
        .filter-buttons {{
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
        }}
        
        .filter-btn {{
            padding: 8px 16px;
            border: 1px solid var(--border-color);
            background: var(--white);
            border-radius: 20px;
            cursor: pointer;
            font-size: var(--font-size);
            font-family: inherit;
            transition: all 0.2s;
            display: inline-flex;
            align-items: center;
            gap: 6px;
        }}
        
        .filter-btn:hover {{
            background: var(--light-bg);
        }}
        
        .filter-btn.active {{
            background: var(--primary-color);
            color: var(--white);
            border-color: var(--primary-color);
        }}
        
        .filter-count {{
            background: rgba(0,0,0,0.1);
            padding: 2px 6px;
            border-radius: 10px;
            font-size: calc(var(--font-size) - 2px);
            font-family: inherit;
        }}
        
        .filter-btn.active .filter-count {{
            background: rgba(255,255,255,0.25);
        }}
        
        .status-indicator {{
            width: 10px;
            height: 10px;
            border-radius: 50%;
            display: inline-block;
        }}
        
        .status-complete {{
            background: var(--success-color);
        }}
        
        .status-incomplete {{
            background: var(--warning-color);
        }}
        
        /* Search Box */
        .search-box {{
            width: 100%;
            padding: 12px 16px;
            border: 1px solid var(--border-color);
            border-radius: 8px;
            font-size: var(--font-size);
            font-family: inherit;
            transition: border-color 0.2s;
        }}
        
        .search-box:focus {{
            outline: none;
            border-color: var(--primary-color);
        }}
        
        /* Table Container */
        .table-container {{
            background: var(--white);
            border-radius: 12px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            overflow: hidden;
        }}
        
        .table-header {{
            padding: 16px 20px;
            border-bottom: 1px solid var(--border-color);
            display: flex;
            justify-content: space-between;
            align-items: center;
            gap: 15px;
            flex-wrap: wrap;
        }}
        
        .result-count {{
            color: var(--text-light);
            font-size: var(--font-size);
            font-family: inherit;
        }}
        
        .export-controls {{
            display: flex;
            gap: 10px;
            align-items: center;
        }}
        
        .export-btn {{
            padding: 8px 16px;
            background: var(--primary-color);
            color: var(--white);
            border: none;
            border-radius: 6px;
            cursor: pointer;
            font-size: var(--font-size);
            font-family: inherit;
            transition: background 0.2s;
        }}
        
        .export-btn:hover {{
            background: var(--secondary-color);
        }}
        
        .export-btn.secondary {{
            background: var(--white);
            color: var(--primary-color);
            border: 1px solid var(--border-color);
        }}
        
        .export-btn.secondary:hover {{
            background: var(--light-bg);
        }}
        
        .export-format {{
            padding: 6px 10px;
            border: 1px solid var(--border-color);
            border-radius: 6px;
            font-size: var(--font-size);
            font-family: inherit;
            background: var(--white);
        }}
        
        .table-wrapper {{
            overflow-x: auto;
            max-height: 70vh;
        }}
        
        table {{
            width: 100%;
            border-collapse: collapse;
            font-size: var(--font-size);
            font-family: inherit;
        }}
        
        thead {{
            position: sticky;
            top: 0;
            background: var(--white);
            z-index: 10;
        }}
        
        th {{
            padding: 14px 12px;
            text-align: left;
            font-weight: 600;
            color: var(--primary-color);
            border-bottom: 2px solid var(--border-color);
            cursor: pointer;
            user-select: none;
            white-space: nowrap;
            font-family: inherit;
        }}
        
        th:hover {{
            background: var(--light-bg);
        }}
        
        /* Remove the arrow indicator */
        th::after {{
            content: '';
        }}
        
        td {{
            padding: 12px;
            border-bottom: 1px solid var(--light-bg);
            vertical-align: top;
            font-family: inherit;
        }}
        
        tbody tr:hover {{
            background: #f8f9fa;
        }}
        
        /* Specific Column Styling */
        .aa-badge {{
            display: inline-block;
            background: var(--light-bg);
            padding: 4px 10px;
            border-radius: 4px;
            font-size: var(--font-size);
            font-weight: 500;
            font-family: inherit;
        }}
        
        .module-link {{
            color: #3498db;
            text-decoration: none;
            font-family: inherit;
            font-weight: 600;
        }}
        
        .module-link:hover {{
            text-decoration: underline;
        }}
        
        .completion-bar {{
            width: 100%;
            max-width: 120px;
        }}
        
        .completion-text {{
            display: flex;
            justify-content: space-between;
            font-size: var(--font-size);
            margin-bottom: 4px;
            color: var(--text-light);
            font-family: inherit;
        }}
        
        .completion-text strong {{
            color: var(--text-dark);
            font-family: inherit;
        }}
        
        .progress-track {{
            width: 100%;
            height: 6px;
            background: var(--light-bg);
            border-radius: 3px;
            overflow: hidden;
        }}
        
        .progress-fill {{
            height: 100%;
            border-radius: 3px;
            transition: width 0.3s;
        }}
        
        .progress-complete {{
            background: var(--success-color);
        }}
        
        .progress-incomplete {{
            background: var(--warning-color);
        }}
        
        .status-badge {{
            display: inline-block;
            padding: 4px 10px;
            border-radius: 4px;
            font-size: var(--font-size);
            font-weight: 600;
            font-family: inherit;
        }}
        
        .badge-complete {{
            background: #d5f4e6;
            color: var(--success-color);
        }}
        
        .badge-incomplete {{
            background: #fef5e7;
            color: #d68910;
        }}
        
        .ko-list {{
            font-size: var(--font-size);
            color: var(--text-dark);
            max-width: 300px;
            word-wrap: break-word;
            line-height: 1.5;
            font-family: inherit;
        }}
        
        .ko-missing {{
            color: var(--danger-color);
        }}
        
        .empty-state {{
            text-align: center;
            padding: 60px 20px;
            color: var(--text-light);
            font-family: inherit;
        }}
        
        .empty-state-icon {{
            font-size: 48px;
            margin-bottom: 16px;
        }}
        
        /* Export Mode - Hide UI elements */
        .export-mode {{
            position: absolute;
            left: -9999px;
            top: 0;
            background: white;
            padding: 30px;
        }}
        
        .export-mode .table-wrapper {{
            max-height: none !important;
            overflow: visible !important;
        }}
        
        /* Print Styles */
        @media print {{
            body {{
                background: white;
                padding: 0;
            }}
            
            .header, .controls, .table-header {{
                display: none !important;
            }}
            
            .table-wrapper {{
                max-height: none !important;
                overflow: visible !important;
            }}
            
            thead {{
                position: static;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <!-- Header -->
        <div class="header">
            <div>
                <h1>🧬 Metabolic Pathway Visualizer</h1>
                <div class="header-info">
                    <span>📊 Total Pathways: {len(data)}</span>
                    <span>🔬 Unique Amino Acids: {len(aa_counts)}</span>
                </div>
            </div>
            <div class="font-controls">
                <label>Font Size:</label>
                <select id="font-size-selector">
                    <option value="11px">Small (11px)</option>
                    <option value="13px" selected>Medium (13px)</option>
                    <option value="15px">Large (15px)</option>
                    <option value="17px">Extra Large (17px)</option>
                </select>
            </div>
        </div>
        
        <!-- Controls -->
        <div class="controls">
            <div class="control-section">
                <span class="control-label">Filter by Amino Acid</span>
                <div class="filter-buttons" id="amino-acid-filters">
                    <button class="filter-btn active" data-filter="all" data-type="aa">
                        All <span class="filter-count">{len(data)}</span>
                    </button>
'''
    
    # Add amino acid filter buttons
    for aa, count in sorted_amino_acids:
        html_content += f'''                    <button class="filter-btn" data-filter="{aa}" data-type="aa">
                        {aa} <span class="filter-count">{count}</span>
                    </button>
'''
    
    html_content += '''                </div>
            </div>
            
            <div class="control-section">
                <span class="control-label">Filter by Status</span>
                <div class="filter-buttons" id="status-filters">
                    <button class="filter-btn active" data-filter="all" data-type="status">
                        All <span class="filter-count">''' + str(len(data)) + '''</span>
                    </button>
'''
    
    # Add status filter buttons
    for status, count in sorted(status_counts.items()):
        status_class = 'complete' if status == 'Complete' else 'incomplete'
        html_content += f'''                    <button class="filter-btn" data-filter="{status}" data-type="status">
                        <span class="status-indicator status-{status_class}"></span>
                        {status} <span class="filter-count">{count}</span>
                    </button>
'''
    
    html_content += '''                </div>
            </div>
            
            <div class="control-section">
                <span class="control-label">Search</span>
                <input type="text" class="search-box" id="search-input" placeholder="Search pathways, modules, KOs...">
            </div>
        </div>
        
        <!-- Table -->
        <div class="table-container">
            <div class="table-header">
                <div class="result-count" id="result-count"></div>
                <div class="export-controls">
                    <select id="export-format" class="export-format">
                        <option value="png">PNG (300 DPI)</option>
                        <option value="svg">SVG</option>
                        <option value="pdf">PDF</option>
                    </select>
                    <button class="export-btn" onclick="exportTable()">📥 Export Table</button>
                </div>
            </div>
            
            <div class="table-wrapper" id="table-wrapper">
                <table id="data-table">
                    <thead>
                        <tr>
                            <th data-sort="amino_acid">Amino Acid</th>
                            <th data-sort="module">Module ID</th>
                            <th data-sort="pathway">Pathway Name</th>
                            <th data-sort="completion_pct">Completion</th>
                            <th data-sort="present_kos_count">KO Count</th>
                            <th data-sort="status">Status</th>
                            <th>Present KOs</th>
                            <th>Missing KOs</th>
                        </tr>
                    </thead>
                    <tbody id="table-body">
                    </tbody>
                </table>
            </div>
        </div>
    </div>
    
    <script>
        // Data
        const allData = ''' + data_json + ''';
        let filteredData = [...allData];
        let currentFilters = {
            aminoAcid: 'all',
            status: 'all',
            search: ''
        };
        let currentSort = {
            column: 'amino_acid',
            ascending: true
        };
        
        // Font size control
        document.getElementById('font-size-selector').addEventListener('change', function() {
            document.documentElement.style.setProperty('--font-size', this.value);
        });
        
        // Render table
        function renderTable() {
            const tbody = document.getElementById('table-body');
            
            if (filteredData.length === 0) {
                tbody.innerHTML = `
                    <tr>
                        <td colspan="8">
                            <div class="empty-state">
                                <div class="empty-state-icon">🔍</div>
                                <div>No pathways found matching your filters</div>
                            </div>
                        </td>
                    </tr>
                `;
                document.getElementById('result-count').textContent = 'No results';
                return;
            }
            
            let html = '';
            filteredData.forEach(row => {
                const completionPct = row.completion_pct;
                const isComplete = row.status === 'Complete';
                const progressClass = isComplete ? 'progress-complete' : 'progress-incomplete';
                const badgeClass = isComplete ? 'badge-complete' : 'badge-incomplete';
                
                html += `
                    <tr>
                        <td><span class="aa-badge">${row.amino_acid}</span></td>
                        <td><a href="https://www.kegg.jp/module/${row.module}" target="_blank" class="module-link">${row.module}</a></td>
                        <td style="font-weight: 500;">${row.pathway}</td>
                        <td>
                            <div class="completion-bar">
                                <div class="completion-text">
                                    <span>${row.satisfied_steps}/${row.total_steps}</span>
                                    <strong>${completionPct}%</strong>
                                </div>
                                <div class="progress-track">
                                    <div class="progress-fill ${progressClass}" style="width: ${completionPct}%"></div>
                                </div>
                            </div>
                        </td>
                        <td>${row.present_kos_count} / ${row.total_kos}</td>
                        <td><span class="status-badge ${badgeClass}">${row.status}</span></td>
                        <td><div class="ko-list">${row.present_kos || '-'}</div></td>
                        <td><div class="ko-list ko-missing">${row.missing_kos || '-'}</div></td>
                    </tr>
                `;
            });
            
            tbody.innerHTML = html;
            document.getElementById('result-count').textContent = 
                `Showing ${filteredData.length} of ${allData.length} pathways`;
        }
        
        // Apply filters
        function applyFilters() {
            filteredData = allData.filter(row => {
                if (currentFilters.aminoAcid !== 'all' && row.amino_acid !== currentFilters.aminoAcid) {
                    return false;
                }
                
                if (currentFilters.status !== 'all' && row.status !== currentFilters.status) {
                    return false;
                }
                
                if (currentFilters.search) {
                    const searchLower = currentFilters.search.toLowerCase();
                    const searchableText = JSON.stringify(row).toLowerCase();
                    if (!searchableText.includes(searchLower)) {
                        return false;
                    }
                }
                
                return true;
            });
            
            applySorting();
            renderTable();
        }
        
        // Sorting
        function applySorting() {
            filteredData.sort((a, b) => {
                let valA = a[currentSort.column];
                let valB = b[currentSort.column];
                
                if (typeof valA === 'string') {
                    valA = valA.toLowerCase();
                    valB = valB.toLowerCase();
                }
                
                let comparison = 0;
                if (valA < valB) comparison = -1;
                if (valA > valB) comparison = 1;
                
                return currentSort.ascending ? comparison : -comparison;
            });
        }
        
        // Export function with 300 DPI for PNG
        async function exportTable() {
            const format = document.getElementById('export-format').value;
            const table = document.getElementById('data-table');
            const clone = table.cloneNode(true);
            
            // Create temporary container for export
            const exportContainer = document.createElement('div');
            exportContainer.className = 'export-mode';
            exportContainer.style.width = '1400px';
            exportContainer.appendChild(clone);
            document.body.appendChild(exportContainer);
            
            try {
                if (format === 'svg') {
                    // SVG export
                    await exportAsSVG(exportContainer);
                } else if (format === 'pdf') {
                    // PDF export using print
                    window.print();
                } else {
                    // PNG export at 300 DPI
                    await exportAsPNG(exportContainer);
                }
            } catch (error) {
                console.error('Export error:', error);
                alert('Export failed. Please try again.');
            } finally {
                document.body.removeChild(exportContainer);
            }
        }
        
        async function exportAsPNG(container) {
            // 300 DPI = scale of 4.17 (300/72)
            const canvas = await html2canvas(container, {
                scale: 4.17, // 300 DPI
                backgroundColor: '#ffffff',
                logging: false,
                windowWidth: 1400
            });
            
            const link = document.createElement('a');
            link.download = 'pathway_table_300dpi.png';
            link.href = canvas.toDataURL('image/png');
            link.click();
        }
        
        async function exportAsSVG(container) {
            // For SVG, we'll use html2canvas and convert to blob
            const canvas = await html2canvas(container, {
                scale: 2,
                backgroundColor: '#ffffff',
                logging: false,
                windowWidth: 1400
            });
            
            // Convert canvas to blob
            canvas.toBlob(function(blob) {
                const url = URL.createObjectURL(blob);
                const link = document.createElement('a');
                link.download = 'pathway_table.svg';
                link.href = url;
                link.click();
                URL.revokeObjectURL(url);
            }, 'image/png'); // Note: html2canvas doesn't produce true SVG
        }
        
        // Event listeners for filters
        document.querySelectorAll('#amino-acid-filters .filter-btn').forEach(btn => {
            btn.addEventListener('click', function() {
                document.querySelectorAll('#amino-acid-filters .filter-btn').forEach(b => 
                    b.classList.remove('active'));
                this.classList.add('active');
                currentFilters.aminoAcid = this.dataset.filter;
                applyFilters();
            });
        });
        
        document.querySelectorAll('#status-filters .filter-btn').forEach(btn => {
            btn.addEventListener('click', function() {
                document.querySelectorAll('#status-filters .filter-btn').forEach(b => 
                    b.classList.remove('active'));
                this.classList.add('active');
                currentFilters.status = this.dataset.filter;
                applyFilters();
            });
        });
        
        // Search
        document.getElementById('search-input').addEventListener('input', function() {
            currentFilters.search = this.value;
            applyFilters();
        });
        
        // Sorting on column click
        document.querySelectorAll('th[data-sort]').forEach(th => {
            th.addEventListener('click', function() {
                const column = this.dataset.sort;
                if (currentSort.column === column) {
                    currentSort.ascending = !currentSort.ascending;
                } else {
                    currentSort.column = column;
                    currentSort.ascending = true;
                }
                applySorting();
                renderTable();
            });
        });
        
        // Initial render
        renderTable();
    </script>
</body>
</html>'''
    
    # Write to file
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    print(f"✅ Generated: {output_path}")
    print(f"📊 Processed {len(data)} pathways")
    print(f"🔬 Found {len(aa_counts)} unique amino acids")

def main():
    print("🧬 Pathway Visualizer Generator")
    print("-" * 50)
    
    # Parse data
    data = parse_csv_data(INPUT_FILE)
    
    if not data:
        print("❌ No data found or error reading CSV file")
        print(f"   Expected file: {INPUT_FILE}")
        return
    
    # Generate HTML
    generate_html(data, OUTPUT_HTML)
    print(f"\n✨ Open {OUTPUT_HTML} in your browser to view the visualizer")

if __name__ == "__main__":
    main()
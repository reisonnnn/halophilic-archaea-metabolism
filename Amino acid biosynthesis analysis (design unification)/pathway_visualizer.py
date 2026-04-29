import csv
import json
import sys
import os
from pathlib import Path
from collections import defaultdict
from datetime import datetime

SCRIPT_DIR = Path(__file__).parent.resolve()
INPUT_FILE = SCRIPT_DIR / "amino_acid_module_details.csv"
OUTPUT_HTML = SCRIPT_DIR / "pathway_visualizer.html"


def parse_csv_data(filepath):
    """Parse the amino acid module details CSV"""
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
                completion_pct = round((satisfied_steps / total_steps * 100), 1) if total_steps > 0 else 0

                item = {
                    'amino_acid': row.get('Amino_Acid', '').strip(),
                    'module': row.get('Module_ID', '').strip(),
                    'pathway': row.get('Pathway_Name', '').strip(),
                    'total_steps': total_steps,
                    'satisfied_steps': satisfied_steps,
                    'completion_pct': completion_pct,
                    'total_kos': int(row.get('Total_KOs', 0)),
                    'present_kos_count': int(row.get('Present_KOs_count', 0)),
                    'status': row.get('Status', '').strip(),
                    'present_kos': row.get('Present_KOs', '').strip(),
                    'missing_kos': row.get('Missing_KOs', '').strip()
                }
                data.append(item)
            except (ValueError, KeyError) as e:
                print(f"Warning: Skipping row: {e}")
                continue
    return data


def generate_html(data, output_path):
    """Generate interactive HTML visualizer matching ko_analyzer design"""

    # Compute counts
    aa_counts = defaultdict(int)
    status_counts = defaultdict(int)
    for item in data:
        aa_counts[item['amino_acid']] += 1
        status_counts[item['status']] += 1

    total = len(data)
    unique_aa = len(aa_counts)
    complete_count = status_counts.get('Complete', 0)
    incomplete_count = status_counts.get('Incomplete', 0)
    no_pathway_count = status_counts.get('No Pathway', 0)

    # Build amino acid filter buttons
    aa_btns = ''
    for aa, cnt in sorted(aa_counts.items()):
        aa_btns += f'<button class="btn" data-v="{aa}">{aa}<span class="cnt">{cnt}</span></button>'

    # Status filter buttons with colored indicators
    status_order = ['Complete', 'Incomplete', 'No Pathway']
    status_colors = {'Complete': '#27ae60', 'Incomplete': '#f39c12', 'No Pathway': '#e74c3c'}
    status_btns = ''
    for s in status_order:
        c = status_counts.get(s, 0)
        color = status_colors.get(s, '#999')
        status_btns += f'<button class="btn" data-v="{s}"><span class="sdot" style="background:{color}"></span>{s}<span class="cnt">{c}</span></button>'

    data_json = json.dumps(data)
    date_str = datetime.now().strftime('%Y-%m-%d %H:%M')

    html = f'''<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<title>Amino Acid Biosynthesis Analyzer</title>
<style>
*{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:system-ui,sans-serif;background:#f0f2f5;padding:20px}}
.container{{max-width:1600px;margin:0 auto}}
header{{text-align:center;padding:30px;background:#fff;border-radius:12px;margin-bottom:20px}}
h1{{color:#2c3e50;margin-bottom:10px}}
.sub{{color:#666;font-style:italic;font-size:13px}}
.legend{{display:flex;justify-content:center;gap:20px;margin-top:15px;flex-wrap:wrap}}
.legend span{{display:flex;align-items:center;gap:5px;font-size:13px}}
.ldot{{width:12px;height:12px;border-radius:50%;display:inline-block}}
.stats{{display:flex;justify-content:center;gap:50px;margin-bottom:20px;flex-wrap:wrap}}
.stat{{text-align:center;background:#fff;padding:20px 30px;border-radius:10px;min-width:120px}}
.stat-val{{font-size:36px;font-weight:bold;color:#2c3e50}}
.stat-lbl{{font-size:12px;color:#666;text-transform:uppercase;margin-top:4px}}
.controls{{display:grid;grid-template-columns:1fr 1fr 300px;gap:20px;margin-bottom:20px}}
@media(max-width:1000px){{.controls{{grid-template-columns:1fr}}}}
.panel{{background:#fff;padding:20px;border-radius:10px}}
.panel h3{{font-size:12px;color:#2c3e50;margin-bottom:12px;text-transform:uppercase;letter-spacing:0.5px}}
.filters{{display:flex;flex-wrap:wrap;gap:8px}}
.btn{{padding:6px 12px;border:2px solid #ddd;border-radius:20px;background:#fff;font-size:12px;cursor:pointer;display:inline-flex;align-items:center;gap:4px;transition:all .15s}}
.btn:hover{{border-color:#2c3e50}}
.btn.active{{background:#2c3e50;color:#fff;border-color:#2c3e50}}
.btn .cnt{{background:rgba(0,0,0,0.1);padding:2px 6px;border-radius:10px;margin-left:4px;font-size:10px}}
.btn.active .cnt{{background:rgba(255,255,255,0.2)}}
.sdot{{width:10px;height:10px;border-radius:50%;display:inline-block}}
input[type=text]{{width:100%;padding:10px;border:2px solid #ddd;border-radius:8px;font-size:14px;font-family:inherit}}
input[type=text]:focus{{outline:none;border-color:#2c3e50}}
.tbl-wrap{{background:#fff;border-radius:12px;overflow:hidden;margin-bottom:20px}}
.tbl-head{{display:flex;justify-content:space-between;align-items:center;padding:15px 20px;background:#f8f9fa;border-bottom:1px solid #eee;flex-wrap:wrap;gap:10px}}
.tbl-title{{font-weight:600;color:#2c3e50}}
.tbl-btns{{display:flex;gap:10px;flex-wrap:wrap}}
.abtn{{padding:8px 16px;border:none;border-radius:6px;font-size:13px;cursor:pointer;font-family:inherit;transition:all .15s}}
.abtn-p{{background:#2c3e50;color:#fff}}
.abtn-p:hover{{background:#1a252f}}
.abtn-s{{background:#e9ecef;color:#2c3e50}}
.abtn-s:hover{{background:#dde1e5}}
.tbl-scroll{{max-height:600px;overflow-y:auto}}
table{{width:100%;border-collapse:collapse;font-size:13px}}
thead{{position:sticky;top:0;z-index:2}}
th{{background:#2c3e50;color:#fff;padding:12px 8px;text-align:left;font-size:11px;text-transform:uppercase;cursor:pointer;white-space:nowrap;user-select:none}}
th:hover{{background:#1a252f}}
td{{padding:10px 8px;border-bottom:1px solid #eee;vertical-align:middle}}
tr:hover td{{background:rgba(0,0,0,0.02)}}
.mono{{font-family:monospace;font-size:12px}}
.center{{text-align:center}}
input[type=checkbox]{{width:16px;height:16px;cursor:pointer}}
.mod-link{{color:#3498db;text-decoration:none;font-family:monospace;font-size:12px}}
.mod-link:hover{{text-decoration:underline}}
.progress-wrap{{width:100%;min-width:120px}}
.progress-info{{display:flex;justify-content:space-between;font-size:11px;margin-bottom:3px}}
.progress-track{{width:100%;height:8px;background:#eee;border-radius:4px;overflow:hidden}}
.progress-fill{{height:100%;border-radius:4px;transition:width .3s}}
.fill-complete{{background:#27ae60}}
.fill-incomplete{{background:#f39c12}}
.fill-none{{background:#e74c3c;width:0}}
.status-badge{{padding:3px 10px;border-radius:12px;font-size:11px;font-weight:600;display:inline-block}}
.badge-complete{{background:#d5f5e3;color:#1e8449}}
.badge-incomplete{{background:#fef9e7;color:#b7950b}}
.badge-none{{background:#fadbd8;color:#c0392b}}
.ko-cell{{font-family:monospace;font-size:11px;color:#555;max-width:200px;word-break:break-all}}
.ko-missing{{color:#c0392b}}
.complete-row{{background:#f0faf4}}
.incomplete-row{{background:#fffdf0}}
.none-row{{background:#fdf2f0}}
.export{{background:#fff;padding:20px;border-radius:12px;margin-bottom:20px}}
#preview{{display:none;background:#fff;padding:20px;border-radius:12px;margin-bottom:20px}}
#preview.show{{display:block}}
.prev-opts{{display:flex;gap:20px;margin-bottom:15px;flex-wrap:wrap;align-items:center}}
.prev-opts select,.prev-opts input{{padding:8px;border:2px solid #ddd;border-radius:6px;font-family:inherit}}
#fig{{background:#fff;padding:20px}}
#fig h4{{text-align:center;margin-bottom:5px;color:#2c3e50}}
#fig .fsub{{text-align:center;color:#666;font-style:italic;font-size:11px;margin-bottom:15px}}
#fig table{{width:100%;border-collapse:collapse;font-size:11px}}
#fig th{{background:#2c3e50;color:#fff;padding:8px;font-size:10px;text-transform:uppercase}}
#fig td{{padding:6px 8px;border:1px solid #ddd}}
#fig .grp-hdr{{background:#34495e;color:#fff;font-weight:bold;font-size:11px}}
#fig .lgnd{{margin-top:15px;display:flex;flex-wrap:wrap;gap:15px;font-size:10px;justify-content:center}}
#fig .lgnd span{{display:flex;align-items:center;gap:4px}}
#fig .lbox{{width:12px;height:12px;border:1px solid #999;border-radius:2px}}
#fig .fig-bar{{height:10px;border-radius:3px;display:inline-block;vertical-align:middle}}
.toast{{position:fixed;bottom:30px;right:30px;background:#2c3e50;color:#fff;padding:12px 24px;border-radius:8px;opacity:0;transition:0.3s;z-index:999;pointer-events:none}}
.toast.show{{opacity:1}}
</style>
<script src="https://cdnjs.cloudflare.com/ajax/libs/html2canvas/1.4.1/html2canvas.min.js"></script>
</head>
<body>
<div class="container">
<header>
<h1>Amino Acid Biosynthesis Analyzer</h1>
<p class="sub">KEGG Module Completeness Analysis &middot; Generated {date_str}</p>
<div class="legend">
<span><span class="ldot" style="background:#27ae60"></span><b>Complete</b></span>
<span><span class="ldot" style="background:#f39c12"></span><b>Incomplete</b></span>
<span><span class="ldot" style="background:#e74c3c"></span><b>No Pathway</b></span>
</div>
</header>

<div class="stats">
<div class="stat"><div class="stat-val" id="s-total">{total}</div><div class="stat-lbl">Total Pathways</div></div>
<div class="stat"><div class="stat-val" id="s-show">{total}</div><div class="stat-lbl">Showing</div></div>
<div class="stat"><div class="stat-val" id="s-sel">0</div><div class="stat-lbl">Selected</div></div>
<div class="stat"><div class="stat-val" style="color:#27ae60">{complete_count}</div><div class="stat-lbl">Complete</div></div>
<div class="stat"><div class="stat-val" style="color:#f39c12">{incomplete_count}</div><div class="stat-lbl">Incomplete</div></div>
<div class="stat"><div class="stat-val" style="color:#e74c3c">{no_pathway_count}</div><div class="stat-lbl">No Pathway</div></div>
</div>

<div class="controls">
<div class="panel"><h3>Amino Acid</h3><div class="filters" id="f-aa"><button class="btn active" data-v="all">All<span class="cnt">{total}</span></button>{aa_btns}</div></div>
<div class="panel"><h3>Status</h3><div class="filters" id="f-status"><button class="btn active" data-v="all">All<span class="cnt">{total}</span></button>{status_btns}</div></div>
<div class="panel"><h3>Search</h3><input type="text" id="search" placeholder="Search KO, pathway, module..."></div>
</div>

<div class="tbl-wrap">
<div class="tbl-head">
<span class="tbl-title">Biosynthesis Pathways</span>
<div class="tbl-btns">
<button class="abtn abtn-s" id="btn-sel">Select Visible</button>
<button class="abtn abtn-s" id="btn-desel">Deselect All</button>
<button class="abtn abtn-p" id="btn-fig">Generate Figure</button>
</div>
</div>
<div class="tbl-scroll">
<table>
<thead><tr>
<th style="width:40px"></th>
<th id="th-aa">Amino Acid</th>
<th id="th-mod">Module</th>
<th id="th-path">Pathway</th>
<th id="th-comp">Completion</th>
<th id="th-kos">KOs</th>
<th id="th-stat">Status</th>
<th>Present KOs</th>
<th>Missing KOs</th>
</tr></thead>
<tbody id="tbody"></tbody>
</table>
</div>
</div>

<div id="preview">
<h3 style="margin-bottom:15px">Figure Preview (<span id="p-cnt">0</span> selected)</h3>
<div class="prev-opts">
<select id="p-mode"><option value="aa">Group by Amino Acid</option><option value="status">Group by Status</option><option value="flat">Flat List</option></select>
<select id="p-font"><option value="9">9px</option><option value="10">10px</option><option value="11" selected>11px</option><option value="12">12px</option></select>
<input type="text" id="p-title" value="Amino Acid Biosynthesis Pathways" style="width:300px">
</div>
<div id="prev-area" style="border:1px solid #ddd;margin-bottom:15px;overflow-x:auto"></div>
<div style="display:flex;gap:10px;flex-wrap:wrap">
<button class="abtn abtn-p" id="btn-png">Save PNG (300 DPI)</button>
<button class="abtn abtn-p" id="btn-tiff">Save TIFF (300 DPI)</button>
<button class="abtn abtn-p" id="btn-svg">Save SVG</button>
<button class="abtn abtn-s" id="btn-csv">Selected CSV</button>
<button class="abtn abtn-s" id="btn-close">Close</button>
</div>
</div>

<div class="export">
<h3 style="margin-bottom:15px">Quick Export</h3>
<div style="display:flex;gap:10px;flex-wrap:wrap">
<button class="abtn abtn-s" id="btn-all">All (CSV)</button>
<button class="abtn abtn-s" id="btn-filt">Filtered (CSV)</button>
<button class="abtn abtn-s" id="btn-json">All (JSON)</button>
</div>
</div>
</div>

<div class="toast" id="toast"></div>

<script>
var DATA = {data_json};
var ROWCLS = {{"Complete":"complete-row","Incomplete":"incomplete-row","No Pathway":"none-row"}};
var STATBG = {{"Complete":"#d5f5e3","Incomplete":"#fef9e7","No Pathway":"#fadbd8"}};
var filt = DATA.slice();
var sel = {{}};
var fAA = "all", fST = "all", fS = "", sCol = "amino_acid", sAsc = true;

function render() {{
    var tb = document.getElementById("tbody");
    var h = "";
    for (var i = 0; i < filt.length; i++) {{
        var r = filt[i];
        var cls = ROWCLS[r.status] || "";
        var chk = sel[r.amino_acid + "|" + r.module] ? " checked" : "";
        var fillCls = r.status === "Complete" ? "fill-complete" : r.status === "Incomplete" ? "fill-incomplete" : "fill-none";
        var badgeCls = r.status === "Complete" ? "badge-complete" : r.status === "Incomplete" ? "badge-incomplete" : "badge-none";
        var modLink = r.module.startsWith("M") ? "<a class='mod-link' href='https://www.kegg.jp/module/" + r.module + "' target='_blank'>" + r.module + "</a>" : "<span class='mono'>" + r.module + "</span>";
        h += "<tr class='" + cls + "'>";
        h += "<td><input type='checkbox'" + chk + " data-key='" + r.amino_acid + "|" + r.module + "'></td>";
        h += "<td style='font-weight:600'>" + r.amino_acid + "</td>";
        h += "<td>" + modLink + "</td>";
        h += "<td>" + r.pathway + "</td>";
        h += "<td><div class='progress-wrap'><div class='progress-info'><span>" + r.satisfied_steps + "/" + r.total_steps + "</span><strong>" + r.completion_pct + "%</strong></div><div class='progress-track'><div class='progress-fill " + fillCls + "' style='width:" + r.completion_pct + "%'></div></div></div></td>";
        h += "<td class='center'>" + r.present_kos_count + " / " + r.total_kos + "</td>";
        h += "<td><span class='status-badge " + badgeCls + "'>" + r.status + "</span></td>";
        h += "<td class='ko-cell'>" + (r.present_kos || "-") + "</td>";
        h += "<td class='ko-cell ko-missing'>" + (r.missing_kos || "-") + "</td>";
        h += "</tr>";
    }}
    tb.innerHTML = h || "<tr><td colspan='9' style='text-align:center;padding:40px;color:#999'>No results</td></tr>";
    tb.querySelectorAll("input[type=checkbox]").forEach(function(cb) {{
        cb.onchange = function() {{
            var key = cb.getAttribute("data-key");
            if (cb.checked) sel[key] = true; else delete sel[key];
            document.getElementById("s-sel").textContent = Object.keys(sel).length;
        }};
    }});
}}

function apply() {{
    filt = DATA.filter(function(r) {{
        if (fAA !== "all" && r.amino_acid !== fAA) return false;
        if (fST !== "all" && r.status !== fST) return false;
        if (fS) {{
            var t = (r.amino_acid + " " + r.module + " " + r.pathway + " " + r.present_kos + " " + r.missing_kos).toLowerCase();
            if (t.indexOf(fS) < 0) return false;
        }}
        return true;
    }});
    filt.sort(function(a, b) {{
        var va = a[sCol], vb = b[sCol];
        if (va === undefined) va = "";
        if (vb === undefined) vb = "";
        if (typeof va === "string") {{ va = va.toLowerCase(); vb = (vb+"").toLowerCase(); }}
        if (va < vb) return sAsc ? -1 : 1;
        if (va > vb) return sAsc ? 1 : -1;
        return 0;
    }});
    render();
    document.getElementById("s-show").textContent = filt.length;
}}

function getSel() {{ return DATA.filter(function(r) {{ return sel[r.amino_acid + "|" + r.module]; }}); }}

function updPrev() {{
    var mode = document.getElementById("p-mode").value;
    var fs = document.getElementById("p-font").value;
    var title = document.getElementById("p-title").value;
    var data = getSel();
    var h = "<div id='fig' style='font-size:" + fs + "px'><h4>" + title + "</h4><p class='fsub'>KEGG Module Completeness | ● Complete ◐ Incomplete ○ No Pathway</p>";

    if (mode === "aa") {{
        var g = {{}};
        data.forEach(function(r) {{ if (!g[r.amino_acid]) g[r.amino_acid] = []; g[r.amino_acid].push(r); }});
        h += "<table><thead><tr><th>Module</th><th>Pathway</th><th>Steps</th><th>Completion</th><th>Status</th><th>Present KOs</th><th>Missing KOs</th></tr></thead><tbody>";
        Object.keys(g).sort().forEach(function(aa) {{
            h += "<tr><td colspan='7' class='grp-hdr'>" + aa + "</td></tr>";
            g[aa].forEach(function(r) {{
                var bg = STATBG[r.status] || "#fff";
                var barW = Math.max(r.completion_pct, 0);
                var barColor = r.status === "Complete" ? "#27ae60" : r.status === "Incomplete" ? "#f39c12" : "#e74c3c";
                h += "<tr style='background:" + bg + "'><td style='font-family:monospace'>" + r.module + "</td><td>" + r.pathway + "</td><td style='text-align:center'>" + r.satisfied_steps + "/" + r.total_steps + "</td><td style='width:120px'><div style='display:flex;align-items:center;gap:6px'><div style='flex:1;height:10px;background:#eee;border-radius:3px;overflow:hidden'><div class='fig-bar' style='width:" + barW + "%;background:" + barColor + ";height:100%;display:block'></div></div><span style='font-size:10px;min-width:32px;text-align:right'>" + r.completion_pct + "%</span></div></td><td style='text-align:center'>" + r.status + "</td><td style='font-family:monospace;font-size:9px'>" + (r.present_kos || "-") + "</td><td style='font-family:monospace;font-size:9px;color:#c0392b'>" + (r.missing_kos || "-") + "</td></tr>";
            }});
        }});
        h += "</tbody></table>";
    }} else if (mode === "status") {{
        var ord = ["Complete","Incomplete","No Pathway"];
        var g = {{}};
        data.forEach(function(r) {{ if (!g[r.status]) g[r.status] = []; g[r.status].push(r); }});
        h += "<table><thead><tr><th>Status</th><th>Amino Acid</th><th>Module</th><th>Pathway</th><th>Steps</th><th>Completion</th></tr></thead><tbody>";
        ord.forEach(function(s) {{
            if (!g[s]) return;
            var first = true;
            g[s].forEach(function(r) {{
                var bg = STATBG[r.status] || "#fff";
                h += "<tr style='background:" + bg + "'><td style='font-weight:" + (first?"bold":"normal") + "'>" + (first?s:"") + "</td><td style='font-weight:600'>" + r.amino_acid + "</td><td style='font-family:monospace'>" + r.module + "</td><td>" + r.pathway + "</td><td style='text-align:center'>" + r.satisfied_steps + "/" + r.total_steps + "</td><td style='text-align:center'>" + r.completion_pct + "%</td></tr>";
                first = false;
            }});
        }});
        h += "</tbody></table>";
    }} else {{
        h += "<table><thead><tr><th>Amino Acid</th><th>Module</th><th>Pathway</th><th>Steps</th><th>Completion</th><th>Status</th><th>Present KOs</th><th>Missing KOs</th></tr></thead><tbody>";
        data.forEach(function(r) {{
            var bg = STATBG[r.status] || "#fff";
            h += "<tr style='background:" + bg + "'><td style='font-weight:600'>" + r.amino_acid + "</td><td style='font-family:monospace'>" + r.module + "</td><td>" + r.pathway + "</td><td style='text-align:center'>" + r.satisfied_steps + "/" + r.total_steps + "</td><td style='text-align:center'>" + r.completion_pct + "%</td><td>" + r.status + "</td><td style='font-family:monospace;font-size:9px'>" + (r.present_kos || "-") + "</td><td style='font-family:monospace;font-size:9px;color:#c0392b'>" + (r.missing_kos || "-") + "</td></tr>";
        }});
        h += "</tbody></table>";
    }}
    h += "<div class='lgnd'><span><span class='lbox' style='background:#d5f5e3'></span>Complete</span><span><span class='lbox' style='background:#fef9e7'></span>Incomplete</span><span><span class='lbox' style='background:#fadbd8'></span>No Pathway</span></div></div>";
    document.getElementById("prev-area").innerHTML = h;
}}

function toast(m) {{ var t = document.getElementById("toast"); t.textContent = m; t.classList.add("show"); setTimeout(function() {{ t.classList.remove("show"); }}, 2500); }}
function dl(c, f, m) {{ var b = new Blob([c], {{type:m}}); var a = document.createElement("a"); a.href = URL.createObjectURL(b); a.download = f; a.click(); }}

// Filter: Amino Acid
document.querySelectorAll("#f-aa .btn").forEach(function(b) {{
    b.onclick = function() {{
        document.querySelectorAll("#f-aa .btn").forEach(function(x){{x.classList.remove("active")}});
        b.classList.add("active"); fAA = b.getAttribute("data-v"); apply();
    }};
}});

// Filter: Status
document.querySelectorAll("#f-status .btn").forEach(function(b) {{
    b.onclick = function() {{
        document.querySelectorAll("#f-status .btn").forEach(function(x){{x.classList.remove("active")}});
        b.classList.add("active"); fST = b.getAttribute("data-v"); apply();
    }};
}});

// Search
document.getElementById("search").oninput = function(e) {{ fS = e.target.value.toLowerCase(); apply(); }};

// Column sorting
var sortMap = {{"th-aa":"amino_acid","th-mod":"module","th-path":"pathway","th-comp":"completion_pct","th-kos":"present_kos_count","th-stat":"status"}};
Object.keys(sortMap).forEach(function(id) {{
    document.getElementById(id).onclick = function() {{
        var col = sortMap[id];
        if (sCol === col) sAsc = !sAsc; else {{ sCol = col; sAsc = true; }}
        apply();
    }};
}});

// Select / Deselect
document.getElementById("btn-sel").onclick = function() {{ filt.forEach(function(r) {{ sel[r.amino_acid + "|" + r.module] = true; }}); document.getElementById("s-sel").textContent = Object.keys(sel).length; render(); toast("Selected " + filt.length); }};
document.getElementById("btn-desel").onclick = function() {{ sel = {{}}; document.getElementById("s-sel").textContent = 0; render(); toast("Cleared"); }};

// Figure preview
document.getElementById("btn-fig").onclick = function() {{ if (Object.keys(sel).length === 0) {{ toast("Select pathways first"); return; }} document.getElementById("preview").classList.add("show"); document.getElementById("p-cnt").textContent = Object.keys(sel).length; updPrev(); }};
document.getElementById("btn-close").onclick = function() {{ document.getElementById("preview").classList.remove("show"); }};
document.getElementById("p-mode").onchange = updPrev;
document.getElementById("p-font").onchange = updPrev;
document.getElementById("p-title").oninput = updPrev;

// Export: PNG 300 DPI
document.getElementById("btn-png").onclick = function() {{ var el = document.getElementById("fig"); if (!el) {{ toast("Preview first"); return; }} toast("Generating PNG..."); html2canvas(el, {{scale:4.17, backgroundColor:"#fff"}}).then(function(c) {{ var a = document.createElement("a"); a.download = "amino_acid_biosynthesis_300dpi.png"; a.href = c.toDataURL("image/png"); a.click(); toast("PNG saved!"); }}); }};

// Export: SVG
document.getElementById("btn-svg").onclick = function() {{ var el = document.getElementById("fig"); if (!el) {{ toast("Preview first"); return; }} var svg = "<svg xmlns='http://www.w3.org/2000/svg' width='" + el.offsetWidth + "' height='" + el.offsetHeight + "'><foreignObject width='100%' height='100%'><div xmlns='http://www.w3.org/1999/xhtml'>" + el.outerHTML + "</div></foreignObject></svg>"; var b = new Blob([svg], {{type:"image/svg+xml"}}); var a = document.createElement("a"); a.download = "amino_acid_biosynthesis.svg"; a.href = URL.createObjectURL(b); a.click(); toast("SVG saved!"); }};

// Export: TIFF 300 DPI
document.getElementById("btn-tiff").onclick = function() {{ var el = document.getElementById("fig"); if (!el) {{ toast("Preview first"); return; }} toast("Generating TIFF..."); html2canvas(el, {{scale:3, backgroundColor:"#fff"}}).then(function(c) {{ var w = c.width, h = c.height; var ctx = c.getContext("2d"); var img = ctx.getImageData(0, 0, w, h); var px = img.data; var stripSize = w * h * 3; var ifdOff = 8; var nTags = 12; var ifdSize = 2 + nTags * 12 + 4; var bitsOff = ifdOff + ifdSize; var dataOff = bitsOff + 6; var resOff = dataOff + stripSize; var fileSize = resOff + 16; var buf = new ArrayBuffer(fileSize); var u8 = new Uint8Array(buf); var dv = new DataView(buf); dv.setUint16(0, 0x4949, true); dv.setUint16(2, 42, true); dv.setUint32(4, ifdOff, true); var off = ifdOff; dv.setUint16(off, nTags, true); off += 2; function tag(t, tp, cnt, v) {{ dv.setUint16(off, t, true); dv.setUint16(off+2, tp, true); dv.setUint32(off+4, cnt, true); dv.setUint32(off+8, v, true); off += 12; }} tag(256, 4, 1, w); tag(257, 4, 1, h); tag(258, 3, 3, bitsOff); tag(259, 3, 1, 1); tag(262, 3, 1, 2); tag(273, 4, 1, dataOff); tag(277, 3, 1, 3); tag(278, 4, 1, h); tag(279, 4, 1, stripSize); tag(282, 5, 1, resOff); tag(283, 5, 1, resOff + 8); tag(296, 3, 1, 2); dv.setUint32(off, 0, true); dv.setUint16(bitsOff, 8, true); dv.setUint16(bitsOff+2, 8, true); dv.setUint16(bitsOff+4, 8, true); var p = dataOff; for (var i = 0; i < w * h; i++) {{ u8[p++] = px[i*4]; u8[p++] = px[i*4+1]; u8[p++] = px[i*4+2]; }} dv.setUint32(resOff, 300, true); dv.setUint32(resOff+4, 1, true); dv.setUint32(resOff+8, 300, true); dv.setUint32(resOff+12, 1, true); var blob = new Blob([buf], {{type:"image/tiff"}}); var a = document.createElement("a"); a.download = "amino_acid_biosynthesis_300dpi.tiff"; a.href = URL.createObjectURL(blob); a.click(); toast("TIFF saved!"); }}); }};

// Export: CSV selected
document.getElementById("btn-csv").onclick = function() {{ var d = getSel(); if (!d.length) {{ toast("No selection"); return; }} var c = "Amino_Acid,Module_ID,Pathway_Name,Total_Steps,Satisfied_Steps,Completion_Pct,Total_KOs,Present_KOs_count,Status,Present_KOs,Missing_KOs\\n"; d.forEach(function(r) {{ c += '"'+r.amino_acid+'","'+r.module+'","'+(r.pathway||"").replace(/"/g,'""')+'",'+r.total_steps+','+r.satisfied_steps+','+r.completion_pct+','+r.total_kos+','+r.present_kos_count+',"'+r.status+'","'+(r.present_kos||"").replace(/"/g,'""')+'","'+(r.missing_kos||"").replace(/"/g,'""')+'"\\n'; }}); dl(c, "selected_pathways.csv", "text/csv"); toast("CSV saved"); }};

// Quick export: All CSV
document.getElementById("btn-all").onclick = function() {{ var c = "Amino_Acid,Module_ID,Pathway_Name,Total_Steps,Satisfied_Steps,Completion_Pct,Total_KOs,Present_KOs_count,Status,Present_KOs,Missing_KOs\\n"; DATA.forEach(function(r) {{ c += '"'+r.amino_acid+'","'+r.module+'","'+(r.pathway||"").replace(/"/g,'""')+'",'+r.total_steps+','+r.satisfied_steps+','+r.completion_pct+','+r.total_kos+','+r.present_kos_count+',"'+r.status+'","'+(r.present_kos||"").replace(/"/g,'""')+'","'+(r.missing_kos||"").replace(/"/g,'""')+'"\\n'; }}); dl(c, "all_pathways.csv", "text/csv"); toast("Exported"); }};

// Quick export: Filtered CSV
document.getElementById("btn-filt").onclick = function() {{ var c = "Amino_Acid,Module_ID,Pathway_Name,Total_Steps,Satisfied_Steps,Completion_Pct,Total_KOs,Present_KOs_count,Status,Present_KOs,Missing_KOs\\n"; filt.forEach(function(r) {{ c += '"'+r.amino_acid+'","'+r.module+'","'+(r.pathway||"").replace(/"/g,'""')+'",'+r.total_steps+','+r.satisfied_steps+','+r.completion_pct+','+r.total_kos+','+r.present_kos_count+',"'+r.status+'","'+(r.present_kos||"").replace(/"/g,'""')+'","'+(r.missing_kos||"").replace(/"/g,'""')+'"\\n'; }}); dl(c, "filtered_pathways.csv", "text/csv"); toast("Exported"); }};

// Quick export: JSON
document.getElementById("btn-json").onclick = function() {{ dl(JSON.stringify(DATA,null,2), "all_pathways.json", "application/json"); toast("Exported"); }};

render();
</script>
</body>
</html>'''

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)
    print(f"✅ Generated: {output_path}")
    print(f"📊 Processed {len(data)} pathways across {unique_aa} amino acids")
    print(f"   Complete: {complete_count} | Incomplete: {incomplete_count} | No Pathway: {no_pathway_count}")


def main():
    print("🧬 Amino Acid Biosynthesis Analyzer")
    print("=" * 50)

    data = parse_csv_data(INPUT_FILE)
    if not data:
        print(f"❌ No data found. Expected: {INPUT_FILE}")
        return

    generate_html(data, str(OUTPUT_HTML))
    print(f"\n✨ Open {OUTPUT_HTML} in your browser")


if __name__ == "__main__":
    main()
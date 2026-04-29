import json
import sys
import os
from pathlib import Path
from collections import defaultdict
from datetime import datetime

SCRIPT_DIR = Path(__file__).parent.resolve()
HS_FILE = SCRIPT_DIR / "HS_GhostKOALA.txt"
HM_FILE = SCRIPT_DIR / "HM_GhostKOALA.txt"
HA_FILE = SCRIPT_DIR / "HA_GhostKOALA.txt"
KO_LIST_FILE = SCRIPT_DIR / "ko_list.txt"
KO_BRITE_FILE = SCRIPT_DIR / "ko_brite.txt"
OUTPUT_HTML = SCRIPT_DIR / "ko_analyzer.html"

def parse_ghostkoala(filepath):
    kos = set()
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2 and parts[1]:
                ko = parts[1].strip()
                if ko.startswith('K') and len(ko) == 6:
                    kos.add(ko)
    return kos

def compute_distributions(hs_kos, hm_kos, ha_kos):
    all_kos = hs_kos | hm_kos | ha_kos
    distributions = {}
    for ko in all_kos:
        in_hs, in_hm, in_ha = ko in hs_kos, ko in hm_kos, ko in ha_kos
        if in_hs and in_hm and in_ha: dist = 'Core'
        elif in_hs and in_hm: dist = 'HS-HM'
        elif in_hs and in_ha: dist = 'HS-HA'
        elif in_hm and in_ha: dist = 'HM-HA'
        elif in_hs: dist = 'HS unique'
        elif in_hm: dist = 'HM unique'
        else: dist = 'HA unique'
        distributions[ko] = {'distribution': dist, 'HS': in_hs, 'HM': in_hm, 'HA': in_ha}
    return distributions

def parse_brite_hierarchy(brite_path):
    ko_to_category = {}
    if not os.path.exists(brite_path):
        return ko_to_category
    current_b = "Other"
    with open(brite_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('B  '):
                txt = line[3:].strip()
                if txt and txt[0].isdigit():
                    parts = txt.split(' ', 1)
                    current_b = parts[1] if len(parts) > 1 else txt
                else:
                    current_b = txt if txt else "Other"
            elif line.startswith('D      '):
                parts = line[7:].strip().split()
                if parts and parts[0].startswith('K') and len(parts[0]) == 6:
                    ko = parts[0]
                    if ko not in ko_to_category:
                        ko_to_category[ko] = current_b
    return ko_to_category

def load_kegg_data(ko_list_path, target_kos, brite_categories):
    kegg_data = {ko: {'gene': '', 'definition': '', 'category': brite_categories.get(ko, 'Other')} for ko in target_kos}
    if not os.path.exists(ko_list_path):
        return kegg_data
    with open(ko_list_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            parts = line.split('\t', 1)
            if len(parts) < 2: continue
            ko = parts[0].replace('ko:', '')
            if ko not in kegg_data: continue
            content = parts[1].strip()
            if '; ' in content:
                gene, definition = content.split('; ', 1)
            else:
                gene, definition = '', content
            if '[EC:' in definition:
                definition = definition.split('[EC:')[0].strip()
            kegg_data[ko]['gene'] = gene
            kegg_data[ko]['definition'] = definition
    return kegg_data

def generate_html(ko_data, distributions, kegg_data, output_path):
    combined = []
    for ko in sorted(ko_data):
        d = distributions.get(ko, {})
        k = kegg_data.get(ko, {})
        combined.append({
            'ko': ko,
            'gene': k.get('gene', ''),
            'definition': k.get('definition', ''),
            'category': k.get('category', 'Other'),
            'distribution': d.get('distribution', ''),
            'HS': d.get('HS', False),
            'HM': d.get('HM', False),
            'HA': d.get('HA', False)
        })
    
    dist_counts = defaultdict(int)
    cat_counts = defaultdict(int)
    for item in combined:
        dist_counts[item['distribution']] += 1
        cat_counts[item['category']] += 1
    
    dist_order = ['HS unique', 'HM unique', 'HA unique', 'HS-HM', 'HS-HA', 'HM-HA', 'Core']
    dist_colors = {'HS unique': '#FADBD8', 'HM unique': '#D5F5E3', 'HA unique': '#D6EAF8', 'HS-HM': '#F5EEF8', 'HS-HA': '#FCF3CF', 'HM-HA': '#D1F2EB', 'Core': '#EAECEE'}
    
    dist_btns = ''
    for d in dist_order:
        c = dist_counts.get(d, 0)
        color = dist_colors.get(d, '#eee')
        dist_btns += f'<button class="btn" data-v="{d}"><span class="cbox" style="background:{color}"></span>{d}<span class="cnt">{c}</span></button>'
    
    cat_btns = ''
    for cat, cnt in sorted(cat_counts.items(), key=lambda x: -x[1]):
        cat_btns += f'<button class="btn" data-v="{cat}">{cat}<span class="cnt">{cnt}</span></button>'
    
    data_json = json.dumps(combined)
    total = len(combined)
    date_str = datetime.now().strftime('%Y-%m-%d %H:%M')
    
    html = f'''<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<title>KO Analyzer</title>
<style>
*{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:system-ui,sans-serif;background:#f0f2f5;padding:20px}}
.container{{max-width:1600px;margin:0 auto}}
header{{text-align:center;padding:30px;background:#fff;border-radius:12px;margin-bottom:20px}}
h1{{color:#2c3e50;margin-bottom:10px}}
.sub{{color:#666;font-style:italic}}
.legend{{display:flex;justify-content:center;gap:20px;margin-top:15px;flex-wrap:wrap}}
.legend span{{display:flex;align-items:center;gap:5px;font-size:13px}}
.dot{{width:12px;height:12px;border-radius:50%}}
.stats{{display:flex;justify-content:center;gap:50px;margin-bottom:20px}}
.stat{{text-align:center}}
.stat-val{{font-size:36px;font-weight:bold;color:#2c3e50}}
.stat-lbl{{font-size:12px;color:#666;text-transform:uppercase}}
.controls{{display:grid;grid-template-columns:1fr 1fr 300px;gap:20px;margin-bottom:20px}}
.panel{{background:#fff;padding:20px;border-radius:10px}}
.panel h3{{font-size:12px;color:#2c3e50;margin-bottom:12px;text-transform:uppercase}}
.filters{{display:flex;flex-wrap:wrap;gap:8px}}
.btn{{padding:6px 12px;border:2px solid #ddd;border-radius:20px;background:#fff;font-size:12px;cursor:pointer}}
.btn:hover{{border-color:#2c3e50}}
.btn.active{{background:#2c3e50;color:#fff;border-color:#2c3e50}}
.btn .cnt{{background:rgba(0,0,0,0.1);padding:2px 6px;border-radius:10px;margin-left:4px;font-size:10px}}
.btn.active .cnt{{background:rgba(255,255,255,0.2)}}
.cbox{{width:14px;height:14px;border-radius:3px;display:inline-block;margin-right:4px;border:1px solid rgba(0,0,0,0.1)}}
input[type=text]{{width:100%;padding:10px;border:2px solid #ddd;border-radius:8px;font-size:14px}}
input[type=text]:focus{{outline:none;border-color:#2c3e50}}
.tbl-wrap{{background:#fff;border-radius:12px;overflow:hidden;margin-bottom:20px}}
.tbl-head{{display:flex;justify-content:space-between;align-items:center;padding:15px 20px;background:#f8f9fa;border-bottom:1px solid #eee}}
.tbl-title{{font-weight:600;color:#2c3e50}}
.tbl-btns{{display:flex;gap:10px}}
.abtn{{padding:8px 16px;border:none;border-radius:6px;font-size:13px;cursor:pointer}}
.abtn-p{{background:#2c3e50;color:#fff}}
.abtn-s{{background:#e9ecef;color:#2c3e50}}
.tbl-scroll{{max-height:600px;overflow-y:auto}}
table{{width:100%;border-collapse:collapse;font-size:13px}}
thead{{position:sticky;top:0}}
th{{background:#2c3e50;color:#fff;padding:12px 8px;text-align:left;font-size:11px;text-transform:uppercase;cursor:pointer}}
th:hover{{background:#1a252f}}
td{{padding:10px 8px;border-bottom:1px solid #eee}}
tr:hover td{{background:rgba(0,0,0,0.02)}}
.mono{{font-family:monospace}}
.gene{{font-style:italic}}
.center{{text-align:center}}
.pdot{{width:12px;height:12px;border-radius:50%;display:inline-block}}
input[type=checkbox]{{width:16px;height:16px}}
.hs-unique{{background:#FADBD8}}
.hm-unique{{background:#D5F5E3}}
.ha-unique{{background:#D6EAF8}}
.hs-hm{{background:#F5EEF8}}
.hs-ha{{background:#FCF3CF}}
.hm-ha{{background:#D1F2EB}}
.core{{background:#EAECEE}}
.export{{background:#fff;padding:20px;border-radius:12px;margin-bottom:20px}}
#preview{{display:none;background:#fff;padding:20px;border-radius:12px;margin-bottom:20px}}
#preview.show{{display:block}}
.prev-opts{{display:flex;gap:20px;margin-bottom:15px;flex-wrap:wrap}}
.prev-opts select,.prev-opts input{{padding:8px;border:2px solid #ddd;border-radius:6px}}
#fig{{background:#fff;padding:20px}}
#fig h4{{text-align:center;margin-bottom:5px}}
#fig .fsub{{text-align:center;color:#666;font-style:italic;font-size:11px;margin-bottom:15px}}
#fig table{{width:100%;border-collapse:collapse;font-size:11px}}
#fig th{{background:#2c3e50;color:#fff;padding:8px;font-size:10px}}
#fig td{{padding:6px;border:1px solid #ccc}}
#fig .cat-hdr{{background:#34495e;color:#fff;font-weight:bold}}
#fig .lgnd{{margin-top:15px;display:flex;flex-wrap:wrap;gap:10px;font-size:10px}}
#fig .lgnd span{{display:flex;align-items:center;gap:4px}}
#fig .lbox{{width:12px;height:12px;border:1px solid #999}}
.toast{{position:fixed;bottom:30px;right:30px;background:#2c3e50;color:#fff;padding:12px 24px;border-radius:8px;opacity:0;transition:0.3s}}
.toast.show{{opacity:1}}
</style>
<script src="https://cdnjs.cloudflare.com/ajax/libs/html2canvas/1.4.1/html2canvas.min.js"></script>
</head>
<body>
<div class="container">
<header>
<h1>KO Distribution Analyzer</h1>

<div class="legend">
<span><span class="dot" style="background:#E74C3C"></span><b>HS</b>: Halobacterium salinarum</span>
<span><span class="dot" style="background:#27AE60"></span><b>HM</b>: Halomicrobium sp.</span>
<span><span class="dot" style="background:#3498DB"></span><b>HA</b>: Haloarcula sp.</span>
</div>
</header>
<div class="stats">
<div class="stat"><div class="stat-val" id="s-total">{total}</div><div class="stat-lbl">Total</div></div>
<div class="stat"><div class="stat-val" id="s-show">{total}</div><div class="stat-lbl">Showing</div></div>
<div class="stat"><div class="stat-val" id="s-sel">0</div><div class="stat-lbl">Selected</div></div>
</div>
<div class="controls">
<div class="panel"><h3>Distribution</h3><div class="filters" id="f-dist"><button class="btn active" data-v="all">All<span class="cnt">{total}</span></button>{dist_btns}</div></div>
<div class="panel"><h3>Category</h3><div class="filters" id="f-cat"><button class="btn active" data-v="all">All<span class="cnt">{total}</span></button>{cat_btns}</div></div>
<div class="panel"><h3>Search</h3><input type="text" id="search" placeholder="Search..."></div>
</div>
<div class="tbl-wrap">
<div class="tbl-head"><span class="tbl-title">KEGG Orthologs</span><div class="tbl-btns"><button class="abtn abtn-s" id="btn-sel">Select Visible</button><button class="abtn abtn-s" id="btn-desel">Deselect All</button><button class="abtn abtn-p" id="btn-fig">Generate Figure</button></div></div>
<div class="tbl-scroll"><table><thead><tr><th style="width:40px"></th><th id="th-ko">KO ID</th><th id="th-gene">Gene</th><th id="th-def">Protein/Enzyme</th><th id="th-cat">Category</th><th id="th-dist">Distribution</th><th class="center">HS</th><th class="center">HM</th><th class="center">HA</th></tr></thead><tbody id="tbody"></tbody></table></div>
</div>
<div id="preview">
<h3 style="margin-bottom:15px">Figure Preview (<span id="p-cnt">0</span> selected)</h3>
<div class="prev-opts">
<select id="p-mode"><option value="cat">Group by Category</option><option value="dist">Group by Distribution</option><option value="flat">Flat List</option></select>
<select id="p-font"><option value="9">9px</option><option value="10">10px</option><option value="11" selected>11px</option><option value="12">12px</option></select>
<input type="text" id="p-title" value="Key KEGG Orthologs" style="width:300px">
</div>
<div id="prev-area" style="border:1px solid #ddd;margin-bottom:15px;overflow-x:auto"></div>
<div style="display:flex;gap:10px"><button class="abtn abtn-p" id="btn-png">Save PNG</button><button class="abtn abtn-p" id="btn-tiff">Save TIFF</button><button class="abtn abtn-p" id="btn-svg">Save SVG</button><button class="abtn abtn-s" id="btn-csv">Save CSV</button><button class="abtn abtn-s" id="btn-close">Close</button></div>
</div>
<div class="export"><h3 style="margin-bottom:15px">Quick Export</h3><div style="display:flex;gap:10px"><button class="abtn abtn-s" id="btn-all">All (CSV)</button><button class="abtn abtn-s" id="btn-filt">Filtered (CSV)</button><button class="abtn abtn-s" id="btn-json">All (JSON)</button></div></div>
</div>
<div class="toast" id="toast"></div>
<script>
var DATA = {data_json};
var DCLS = {{"HS unique":"hs-unique","HM unique":"hm-unique","HA unique":"ha-unique","HS-HM":"hs-hm","HS-HA":"hs-ha","HM-HA":"hm-ha","Core":"core"}};
var DBGC = {{"HS unique":"#FADBD8","HM unique":"#D5F5E3","HA unique":"#D6EAF8","HS-HM":"#F5EEF8","HS-HA":"#FCF3CF","HM-HA":"#D1F2EB","Core":"#EAECEE"}};
var filt = DATA.slice();
var sel = {{}};
var fD = "all", fC = "all", fS = "", sCol = "ko", sAsc = true;

function render() {{
    var tb = document.getElementById("tbody");
    var h = "";
    for (var i = 0; i < filt.length; i++) {{
        var r = filt[i];
        var cls = DCLS[r.distribution] || "";
        var chk = sel[r.ko] ? " checked" : "";
        h += "<tr class='" + cls + "'>";
        h += "<td><input type='checkbox'" + chk + " data-ko='" + r.ko + "'></td>";
        h += "<td class='mono'>" + r.ko + "</td>";
        h += "<td class='gene'>" + (r.gene || "-") + "</td>";
        h += "<td>" + (r.definition || "-") + "</td>";
        h += "<td>" + r.category + "</td>";
        h += "<td>" + r.distribution + "</td>";
        h += "<td class='center'>" + (r.HS ? "<span class='pdot' style='background:#E74C3C'></span>" : "") + "</td>";
        h += "<td class='center'>" + (r.HM ? "<span class='pdot' style='background:#27AE60'></span>" : "") + "</td>";
        h += "<td class='center'>" + (r.HA ? "<span class='pdot' style='background:#3498DB'></span>" : "") + "</td>";
        h += "</tr>";
    }}
    tb.innerHTML = h || "<tr><td colspan='9' style='text-align:center;padding:40px'>No results</td></tr>";
    
    tb.querySelectorAll("input[type=checkbox]").forEach(function(cb) {{
        cb.onchange = function() {{
            var ko = cb.getAttribute("data-ko");
            if (cb.checked) sel[ko] = true; else delete sel[ko];
            document.getElementById("s-sel").textContent = Object.keys(sel).length;
        }};
    }});
}}

function apply() {{
    filt = DATA.filter(function(r) {{
        if (fD !== "all" && r.distribution !== fD) return false;
        if (fC !== "all" && r.category !== fC) return false;
        if (fS && (r.ko + " " + r.gene + " " + r.definition).toLowerCase().indexOf(fS) < 0) return false;
        return true;
    }});
    filt.sort(function(a, b) {{
        var va = a[sCol] || "", vb = b[sCol] || "";
        if (typeof va === "string") va = va.toLowerCase();
        if (typeof vb === "string") vb = vb.toLowerCase();
        if (va < vb) return sAsc ? -1 : 1;
        if (va > vb) return sAsc ? 1 : -1;
        return 0;
    }});
    render();
    document.getElementById("s-show").textContent = filt.length;
}}

function getSel() {{ return DATA.filter(function(r) {{ return sel[r.ko]; }}); }}

function updPrev() {{
    var mode = document.getElementById("p-mode").value;
    var fs = document.getElementById("p-font").value;
    var title = document.getElementById("p-title").value;
    var data = getSel();
    var h = "<div id='fig' style='font-size:" + fs + "px'><h4>" + title + "</h4><p class='fsub'>HS: Halobacterium salinarum | HM: Halomicrobium sp. | HA: Haloarcula sp.</p>";
    if (mode === "cat") {{
        var g = {{}};
        data.forEach(function(r) {{ if (!g[r.category]) g[r.category] = []; g[r.category].push(r); }});
        h += "<table><thead><tr><th>KO</th><th>Gene</th><th>Protein/Enzyme</th><th>HS</th><th>HM</th><th>HA</th></tr></thead><tbody>";
        Object.keys(g).sort().forEach(function(c) {{
            h += "<tr><td colspan='6' class='cat-hdr'>" + c + "</td></tr>";
            g[c].forEach(function(r) {{
                h += "<tr style='background:" + DBGC[r.distribution] + "'><td style='font-family:monospace'>" + r.ko + "</td><td style='font-style:italic'>" + (r.gene||"-") + "</td><td>" + (r.definition||"-") + "</td><td style='text-align:center;color:#E74C3C'>" + (r.HS?"●":"") + "</td><td style='text-align:center;color:#27AE60'>" + (r.HM?"●":"") + "</td><td style='text-align:center;color:#3498DB'>" + (r.HA?"●":"") + "</td></tr>";
            }});
        }});
        h += "</tbody></table>";
    }} else if (mode === "dist") {{
        var ord = ["HS unique","HM unique","HA unique","HS-HM","HS-HA","HM-HA","Core"];
        var g = {{}};
        data.forEach(function(r) {{ if (!g[r.distribution]) g[r.distribution] = []; g[r.distribution].push(r); }});
        h += "<table><thead><tr><th>Distribution</th><th>KO</th><th>Gene</th><th>Protein/Enzyme</th><th>Category</th></tr></thead><tbody>";
        ord.forEach(function(d) {{
            if (!g[d]) return;
            var first = true;
            g[d].forEach(function(r) {{
                h += "<tr style='background:" + DBGC[r.distribution] + "'><td style='font-weight:" + (first?"bold":"normal") + "'>" + (first?d:"") + "</td><td style='font-family:monospace'>" + r.ko + "</td><td style='font-style:italic'>" + (r.gene||"-") + "</td><td>" + (r.definition||"-") + "</td><td>" + r.category + "</td></tr>";
                first = false;
            }});
        }});
        h += "</tbody></table>";
    }} else {{
        h += "<table><thead><tr><th>KO</th><th>Gene</th><th>Protein/Enzyme</th><th>Category</th><th>Distribution</th><th>HS</th><th>HM</th><th>HA</th></tr></thead><tbody>";
        data.forEach(function(r) {{
            h += "<tr style='background:" + DBGC[r.distribution] + "'><td style='font-family:monospace'>" + r.ko + "</td><td style='font-style:italic'>" + (r.gene||"-") + "</td><td>" + (r.definition||"-") + "</td><td>" + r.category + "</td><td>" + r.distribution + "</td><td style='text-align:center;color:#E74C3C'>" + (r.HS?"●":"") + "</td><td style='text-align:center;color:#27AE60'>" + (r.HM?"●":"") + "</td><td style='text-align:center;color:#3498DB'>" + (r.HA?"●":"") + "</td></tr>";
        }});
        h += "</tbody></table>";
    }}
    h += "<div class='lgnd'><span><span class='lbox' style='background:#FADBD8'></span>HS unique</span><span><span class='lbox' style='background:#D5F5E3'></span>HM unique</span><span><span class='lbox' style='background:#D6EAF8'></span>HA unique</span><span><span class='lbox' style='background:#F5EEF8'></span>HS-HM</span><span><span class='lbox' style='background:#FCF3CF'></span>HS-HA</span><span><span class='lbox' style='background:#D1F2EB'></span>HM-HA</span><span><span class='lbox' style='background:#EAECEE'></span>Core</span></div></div>";
    document.getElementById("prev-area").innerHTML = h;
}}

function toast(m) {{ var t = document.getElementById("toast"); t.textContent = m; t.classList.add("show"); setTimeout(function() {{ t.classList.remove("show"); }}, 2500); }}

function dl(c, f, m) {{ var b = new Blob([c], {{type:m}}); var a = document.createElement("a"); a.href = URL.createObjectURL(b); a.download = f; a.click(); }}

document.querySelectorAll("#f-dist .btn").forEach(function(b) {{
    b.onclick = function() {{
        document.querySelectorAll("#f-dist .btn").forEach(function(x){{x.classList.remove("active")}});
        b.classList.add("active");
        fD = b.getAttribute("data-v");
        apply();
    }};
}});

document.querySelectorAll("#f-cat .btn").forEach(function(b) {{
    b.onclick = function() {{
        document.querySelectorAll("#f-cat .btn").forEach(function(x){{x.classList.remove("active")}});
        b.classList.add("active");
        fC = b.getAttribute("data-v");
        apply();
    }};
}});

document.getElementById("search").oninput = function(e) {{ fS = e.target.value.toLowerCase(); apply(); }};

document.getElementById("th-ko").onclick = function() {{ if (sCol === "ko") sAsc = !sAsc; else {{ sCol = "ko"; sAsc = true; }} apply(); }};
document.getElementById("th-gene").onclick = function() {{ if (sCol === "gene") sAsc = !sAsc; else {{ sCol = "gene"; sAsc = true; }} apply(); }};
document.getElementById("th-def").onclick = function() {{ if (sCol === "definition") sAsc = !sAsc; else {{ sCol = "definition"; sAsc = true; }} apply(); }};
document.getElementById("th-cat").onclick = function() {{ if (sCol === "category") sAsc = !sAsc; else {{ sCol = "category"; sAsc = true; }} apply(); }};
document.getElementById("th-dist").onclick = function() {{ if (sCol === "distribution") sAsc = !sAsc; else {{ sCol = "distribution"; sAsc = true; }} apply(); }};

document.getElementById("btn-sel").onclick = function() {{ filt.forEach(function(r) {{ sel[r.ko] = true; }}); document.getElementById("s-sel").textContent = Object.keys(sel).length; render(); toast("Selected " + filt.length); }};
document.getElementById("btn-desel").onclick = function() {{ sel = {{}}; document.getElementById("s-sel").textContent = 0; render(); toast("Cleared"); }};
document.getElementById("btn-fig").onclick = function() {{ if (Object.keys(sel).length === 0) {{ toast("Select KOs first"); return; }} document.getElementById("preview").classList.add("show"); document.getElementById("p-cnt").textContent = Object.keys(sel).length; updPrev(); }};
document.getElementById("btn-close").onclick = function() {{ document.getElementById("preview").classList.remove("show"); }};

document.getElementById("p-mode").onchange = updPrev;
document.getElementById("p-font").onchange = updPrev;
document.getElementById("p-title").oninput = updPrev;

document.getElementById("btn-png").onclick = function() {{ var el = document.getElementById("fig"); if (!el) {{ toast("Preview first"); return; }} toast("Generating..."); html2canvas(el, {{scale:3, backgroundColor:"#fff"}}).then(function(c) {{ var a = document.createElement("a"); a.download = "ko_table.png"; a.href = c.toDataURL("image/png"); a.click(); toast("PNG saved!"); }}); }};
document.getElementById("btn-svg").onclick = function() {{ var el = document.getElementById("fig"); if (!el) {{ toast("Preview first"); return; }} var svg = "<svg xmlns='http://www.w3.org/2000/svg' width='" + el.offsetWidth + "' height='" + el.offsetHeight + "'><foreignObject width='100%' height='100%'><div xmlns='http://www.w3.org/1999/xhtml'>" + el.outerHTML + "</div></foreignObject></svg>"; var b = new Blob([svg], {{type:"image/svg+xml"}}); var a = document.createElement("a"); a.download = "ko_table.svg"; a.href = URL.createObjectURL(b); a.click(); toast("SVG saved!"); }};
document.getElementById("btn-tiff").onclick = function() {{ var el = document.getElementById("fig"); if (!el) {{ toast("Preview first"); return; }} toast("Generating TIFF..."); html2canvas(el, {{scale:3, backgroundColor:"#fff"}}).then(function(c) {{ var w = c.width, h = c.height; var ctx = c.getContext("2d"); var img = ctx.getImageData(0, 0, w, h); var px = img.data; var stripSize = w * h * 3; var ifdOff = 8; var nTags = 12; var ifdSize = 2 + nTags * 12 + 4; var bitsOff = ifdOff + ifdSize; var dataOff = bitsOff + 6; var resOff = dataOff + stripSize; var fileSize = resOff + 16; var buf = new ArrayBuffer(fileSize); var u8 = new Uint8Array(buf); var dv = new DataView(buf); dv.setUint16(0, 0x4949, true); dv.setUint16(2, 42, true); dv.setUint32(4, ifdOff, true); var off = ifdOff; dv.setUint16(off, nTags, true); off += 2; function tag(t, tp, cnt, v) {{ dv.setUint16(off, t, true); dv.setUint16(off+2, tp, true); dv.setUint32(off+4, cnt, true); dv.setUint32(off+8, v, true); off += 12; }} tag(256, 4, 1, w); tag(257, 4, 1, h); tag(258, 3, 3, bitsOff); tag(259, 3, 1, 1); tag(262, 3, 1, 2); tag(273, 4, 1, dataOff); tag(277, 3, 1, 3); tag(278, 4, 1, h); tag(279, 4, 1, stripSize); tag(282, 5, 1, resOff); tag(283, 5, 1, resOff + 8); tag(296, 3, 1, 2); dv.setUint32(off, 0, true); dv.setUint16(bitsOff, 8, true); dv.setUint16(bitsOff+2, 8, true); dv.setUint16(bitsOff+4, 8, true); var p = dataOff; for (var i = 0; i < w * h; i++) {{ u8[p++] = px[i*4]; u8[p++] = px[i*4+1]; u8[p++] = px[i*4+2]; }} dv.setUint32(resOff, 300, true); dv.setUint32(resOff+4, 1, true); dv.setUint32(resOff+8, 300, true); dv.setUint32(resOff+12, 1, true); var blob = new Blob([buf], {{type:"image/tiff"}}); var a = document.createElement("a"); a.download = "ko_table.tiff"; a.href = URL.createObjectURL(blob); a.click(); toast("TIFF saved!"); }}); }};
document.getElementById("btn-csv").onclick = function() {{ var d = getSel(); if (!d.length) {{ toast("No selection"); return; }} var c = "KO,Gene,Definition,Category,Distribution,HS,HM,HA\\n"; d.forEach(function(r) {{ c += '"'+r.ko+'","'+r.gene+'","'+(r.definition||"").replace(/"/g,'""')+'","'+r.category+'","'+r.distribution+'",'+r.HS+','+r.HM+','+r.HA+"\\n"; }}); dl(c, "selected.csv", "text/csv"); toast("CSV saved"); }};

document.getElementById("btn-all").onclick = function() {{ var c = "KO,Gene,Definition,Category,Distribution,HS,HM,HA\\n"; DATA.forEach(function(r) {{ c += '"'+r.ko+'","'+r.gene+'","'+(r.definition||"").replace(/"/g,'""')+'","'+r.category+'","'+r.distribution+'",'+r.HS+','+r.HM+','+r.HA+"\\n"; }}); dl(c, "all_kos.csv", "text/csv"); toast("Exported"); }};
document.getElementById("btn-filt").onclick = function() {{ var c = "KO,Gene,Definition,Category,Distribution,HS,HM,HA\\n"; filt.forEach(function(r) {{ c += '"'+r.ko+'","'+r.gene+'","'+(r.definition||"").replace(/"/g,'""')+'","'+r.category+'","'+r.distribution+'",'+r.HS+','+r.HM+','+r.HA+"\\n"; }}); dl(c, "filtered.csv", "text/csv"); toast("Exported"); }};
document.getElementById("btn-json").onclick = function() {{ dl(JSON.stringify(DATA,null,2), "all_kos.json", "application/json"); toast("Exported"); }};

render();
</script>
</body>
</html>'''
    
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)
    print(f"Generated: {output_path}")

def main():
    print("KO Distribution Analyzer")
    print("=" * 40)
    
    for f in [HS_FILE, HM_FILE, HA_FILE]:
        if not f.exists():
            print(f"ERROR: {f.name} not found")
            sys.exit(1)
    
    hs = parse_ghostkoala(HS_FILE)
    hm = parse_ghostkoala(HM_FILE)
    ha = parse_ghostkoala(HA_FILE)
    print(f"HS:{len(hs)} HM:{len(hm)} HA:{len(ha)}")
    
    dist = compute_distributions(hs, hm, ha)
    all_kos = list(dist.keys())
    print(f"Total: {len(all_kos)}")
    
    brite = parse_brite_hierarchy(KO_BRITE_FILE)
    print(f"BRITE mapped: {len(brite)}")
    
    kegg = load_kegg_data(KO_LIST_FILE, all_kos, brite)
    
    generate_html(all_kos, dist, kegg, str(OUTPUT_HTML))
    print("Done! Open: ko_analyzer.html")

if __name__ == '__main__':
    main()
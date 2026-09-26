#!/usr/bin/env python3
"""Render training_data_selection_methods.md as a self-contained HTML artifact page."""
import html, os, re, sys

SRC = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/training_data_selection_methods.md"
OUT = sys.argv[1]
# Sections 8-10 (evidence and alignment, experiment A interim) are written by the
# REVIEW session directly as HTML. They are appended unchanged so a rebuild keeps them.
EXTRA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sections_8-10_REVIEW.html")
lines = open(SRC).read().split("\n")


def inline(t):
    t = html.escape(t, quote=False)
    t = re.sub(r"`([^`]+)`", r"<code>\1</code>", t)
    t = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", t)
    t = re.sub(r"(?<![\w*])\*([^*\n]+?)\*(?![\w*])", r"<em>\1</em>", t)
    return t


def slug(t):
    return re.sub(r"[^a-z0-9]+", "-", t.lower()).strip("-")[:60]


body, toc = [], []
i = 0
list_stack = []  # (indent, tag)


def close_lists(to_indent=-1):
    while list_stack and list_stack[-1][0] > to_indent:
        body.append("</li></%s>" % list_stack.pop()[1])


title = ""
para = []


def flush_para():
    if para:
        body.append("<p>%s</p>" % inline(" ".join(para)))
        para.clear()


while i < len(lines):
    ln = lines[i]
    s = ln.strip()
    if not s:
        flush_para(); close_lists(); i += 1; continue
    m = re.match(r"^(#{1,3}) (.*)", ln)
    if m:
        flush_para(); close_lists()
        lvl, txt = len(m.group(1)), m.group(2)
        if lvl == 1:
            title = txt
        else:
            sid = slug(txt)
            tag = "h2" if lvl == 2 else "h3"
            body.append(f'<{tag} id="{sid}">{inline(txt)}</{tag}>')
            if lvl == 2:
                toc.append((sid, re.sub(r"^\d+\.\s*", "", txt)))
        i += 1; continue
    if s.startswith("|") and i + 1 < len(lines) and re.match(r"^\s*\|[-| :]+\|\s*$", lines[i + 1]):
        flush_para(); close_lists()
        head = [c.strip() for c in s.strip("|").split("|")]
        rows = []
        i += 2
        while i < len(lines) and lines[i].strip().startswith("|"):
            rows.append([c.strip() for c in lines[i].strip().strip("|").split("|")])
            i += 1
        t = ['<div class="tbl"><table><thead><tr>' + "".join(f"<th>{inline(h)}</th>" for h in head) + "</tr></thead><tbody>"]
        for r in rows:
            t.append("<tr>" + "".join(f"<td>{inline(c)}</td>" for c in r) + "</tr>")
        t.append("</tbody></table></div>")
        body.append("".join(t))
        continue
    m = re.match(r"^(\s*)(\d+\.|-) (.*)", ln)
    if m:
        flush_para()
        ind = len(m.group(1)); tag = "ol" if m.group(2)[0].isdigit() else "ul"
        if not list_stack or ind > list_stack[-1][0]:
            body.append(f"<{tag}><li>")
            list_stack.append((ind, tag))
        else:
            close_lists(ind)
            if list_stack and list_stack[-1][0] == ind:
                body.append("</li><li>")
            else:
                body.append(f"<{tag}><li>"); list_stack.append((ind, tag))
        body.append(inline(m.group(3)))
        i += 1; continue
    if list_stack and ln.startswith(" "):
        body.append(" " + inline(s)); i += 1; continue
    close_lists()
    para.append(s)
    i += 1
flush_para(); close_lists()

# short contents labels used by REVIEW on the live page (version 2)
TOC_SHORT = {"8-training-versus-evidence": "Training versus evidence",
             "9-alignment-and-pasa-evidence": "Alignment and PASA evidence",
             "10-gate-calibration": "Calibrating the gate"}
extra_html = open(EXTRA).read() if os.path.exists(EXTRA) else ""
for sid, txt in re.findall(r'<h2 id="([^"]+)">(.*?)</h2>', extra_html):
    short = re.sub(r"^\d+\.\s*", "", re.sub(r"<[^>]+>", "", txt))
    toc.append((sid, TOC_SHORT.get(sid, html.unescape(short))))
toc_html = "".join(f'<li><a href="#{sid}">{html.escape(t)}</a></li>' for sid, t in toc)
page = f"""<title>funannotate Training Selection</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans+Condensed:wght@500;600;700&family=Source+Serif+4:ital,opsz,wght@0,8..60,400;0,8..60,600;1,8..60,400&display=swap">
<style>
:root {{
  --paper: #f6f7f3; --ink: #1c2420; --muted: #5a6660; --accent: #1e6a5a;
  --rule: #d6dcd6; --wash: #eaefea; --code-bg: #e7ece7; --row: #fbfcfa;
  --serif: "Source Serif 4", Georgia, "Times New Roman", serif;
  --sans: "IBM Plex Sans Condensed", "Arial Narrow", Arial, sans-serif;
  --mono: "IBM Plex Mono", ui-monospace, SFMono-Regular, Menlo, monospace;
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    --paper: #111614; --ink: #e2e8e4; --muted: #95a29b; --accent: #5fb39f;
    --rule: #2b3430; --wash: #18201d; --code-bg: #1d2622; --row: #141a17; color-scheme: dark;
  }}
}}
:root[data-theme="dark"] {{
  --paper: #111614; --ink: #e2e8e4; --muted: #95a29b; --accent: #5fb39f;
  --rule: #2b3430; --wash: #18201d; --code-bg: #1d2622; --row: #141a17; color-scheme: dark;
}}
body {{ background: var(--paper); color: var(--ink); font: 17px/1.62 var(--serif); padding-inline: 20px; padding-block: 0 64px; }}
.wrap {{ max-width: 1120px; margin: 0 auto; display: grid; grid-template-columns: minmax(0, 1fr); gap: 32px; }}
@media (min-width: 1000px) {{ .wrap {{ grid-template-columns: 220px minmax(0, 1fr); gap: 48px; }} }}
header.masthead {{ grid-column: 1 / -1; padding-block: 44px 20px; border-bottom: 1px solid var(--rule); }}
.eyebrow {{ font: 600 12px/1.2 var(--sans); letter-spacing: .12em; text-transform: uppercase; color: var(--accent); }}
h1 {{ font: 700 clamp(30px, 4.4vw, 46px)/1.08 var(--sans); margin: 10px 0 12px; text-wrap: balance; letter-spacing: -.01em; }}
.lede {{ max-width: 66ch; color: var(--muted); margin: 0; }}
.meta {{ display: flex; flex-wrap: wrap; gap: 8px 18px; margin-top: 16px; font: 13px/1.4 var(--mono); color: var(--muted); }}
nav.toc {{ align-self: start; }}
@media (min-width: 1000px) {{ nav.toc {{ position: sticky; top: calc(env(safe-area-inset-top, 0px) + 24px); padding-top: 28px; }} }}
nav.toc h2 {{ font: 600 12px var(--sans); letter-spacing: .12em; text-transform: uppercase; color: var(--muted); margin: 0 0 10px; }}
nav.toc ol {{ list-style: none; padding: 0; margin: 0; display: grid; gap: 6px; font: 500 14px/1.3 var(--sans); }}
nav.toc a {{ color: var(--ink); text-decoration: none; border-left: 2px solid var(--rule); padding-left: 10px; display: block; }}
nav.toc a:hover, nav.toc a:focus-visible {{ border-left-color: var(--accent); color: var(--accent); outline: none; }}
main {{ min-width: 0; }}
main > p, main > ul, main > ol {{ max-width: 68ch; }}
main h2 {{ font: 600 26px/1.2 var(--sans); margin: 44px 0 12px; padding-top: 8px; border-top: 1px solid var(--rule); text-wrap: balance; }}
main h3 {{ font: 600 18px/1.3 var(--sans); margin: 32px 0 8px; color: var(--accent); text-wrap: balance; }}
p {{ margin: 0 0 14px; }}
ul, ol {{ margin: 0 0 16px; padding-left: 1.3em; }}
li {{ margin: 4px 0; }}
li > ul, li > ol {{ margin: 6px 0 4px; }}
code {{ font: 0.84em var(--mono); background: var(--code-bg); padding: 1px 5px; border-radius: 3px; word-break: break-word; }}
strong {{ font-weight: 600; }}
.tbl {{ overflow-x: auto; margin: 10px 0 20px; border: 1px solid var(--rule); border-radius: 4px; background: var(--row); }}
table {{ border-collapse: collapse; width: 100%; font: 13.5px/1.45 var(--mono); font-variant-numeric: tabular-nums; }}
th {{ font: 600 12.5px/1.3 var(--sans); letter-spacing: .02em; text-align: left; background: var(--wash); color: var(--ink); padding: 9px 12px; border-bottom: 1px solid var(--rule); vertical-align: bottom; }}
td {{ padding: 7px 12px; border-bottom: 1px solid var(--rule); vertical-align: top; white-space: nowrap; }}
td:first-child {{ font-family: var(--sans); font-size: 14px; white-space: normal; min-width: 9ch; }}
tbody tr:last-child td {{ border-bottom: 0; }}
tbody tr:hover td {{ background: var(--wash); }}
td code {{ background: none; padding: 0; }}
a {{ color: var(--accent); }}
@media (prefers-reduced-motion: reduce) {{ * {{ scroll-behavior: auto !important; }} }}
html {{ scroll-behavior: smooth; }}
</style>
<div class="wrap">
  <header class="masthead">
    <div class="eyebrow">funannotate · methods draft · 2026-09-26</div>
    <h1>{inline(title)}</h1>
    <p class="lede">How funannotate chooses the gene models that train Augustus and snap, which thresholds it applies, how each decision is logged, and how every change was measured against RefSeq on held-out chromosomes of five fungal genomes.</p>
    <div class="meta"><span>funannotate-live 66c75ea</span><span>funannotate-1.9.0-rc.1 image</span><span>5 RefSeq genomes</span><span>tables M1–M8 generated from result files</span><span>§8–10 evidence and alignment: REVIEW session, 2026-09-26</span></div>
  </header>
  <nav class="toc" aria-label="Contents"><h2>Contents</h2><ol>{toc_html}</ol></nav>
  <main>
{chr(10).join(body)}

{extra_html}  </main>
</div>
"""
open(OUT, "w").write(page)
print("wrote", OUT, len(page), "bytes;", len(toc), "sections")

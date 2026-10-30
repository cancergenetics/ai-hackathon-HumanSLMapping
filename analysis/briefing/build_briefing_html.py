"""Render BRIEFING.md to one self-contained HTML page (figures embedded). Run: uv run --with markdown python analysis/briefing/build_briefing_html.py OUT.html"""
import base64, re, sys
from pathlib import Path
import markdown
root = Path(__file__).resolve().parents[2]
md = (root / "BRIEFING.md").read_text()
def img(m):
    alt, path = m.group(1), m.group(2)
    b = base64.b64encode((root / path).read_bytes()).decode()
    return f'![{alt}](data:image/png;base64,{b})'
md = re.sub(r'!\[([^\]]*)\]\(([^)]+\.png)\)', img, md)
md = re.sub(r'!\[([^\]]*)\]\((data:[^)]+)\)\n\*(Figure[^\n]*)\*',
            lambda m: f'<figure><img alt="{m.group(1)}" src="{m.group(2)}"><figcaption>{m.group(3)}</figcaption></figure>', md)
body = markdown.markdown(md, extensions=["tables", "sane_lists"])
# wrap figure + caption: <p><img ...></p>\n<p><em>Figure..</em></p>
body = re.sub(r'<p>(<img [^>]+>)</p>\s*<p><em>(Figure[^<]*(?:<[^/][^>]*>[^<]*</[^>]+>[^<]*)*)</em></p>',
              r'<figure>\1<figcaption>\2</figcaption></figure>', body)
body = re.sub(r'<table>', '<div class="tw"><table>', body); body = body.replace('</table>', '</table></div>')
head = """<title>Choosing Which Genes to Screen</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Source+Serif+4:ital,opsz,wght@0,8..60,400;0,8..60,600;1,8..60,400&family=Source+Sans+3:wght@400;600&family=IBM+Plex+Mono:wght@400&display=swap">
<style>
/* Layout: single journal-style reading column; figures on white plates so the plotted colours stay faithful in both themes. */
:root { --bg:#f5f7f6; --ink:#14211d; --muted:#4f5d59; --line:#d5dedb; --accent:#2a78d6; --plate:#ffffff; --tablehead:#e9efed;
  --display:"Source Serif 4", Georgia, serif; --body:"Source Sans 3", system-ui, -apple-system, "Segoe UI", sans-serif; --mono:"IBM Plex Mono", ui-monospace, Menlo, monospace; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { --bg:#0f1513; --ink:#e6edea; --muted:#a6b3af; --line:#2c3835; --accent:#5b9be9; --plate:#ffffff; --tablehead:#18211f; color-scheme:dark; } }
:root[data-theme="dark"] { --bg:#0f1513; --ink:#e6edea; --muted:#a6b3af; --line:#2c3835; --accent:#5b9be9; --plate:#ffffff; --tablehead:#18211f; color-scheme:dark; }
body { background:var(--bg); color:var(--ink); font-family:var(--body); font-size:17px; line-height:1.62; }
main { max-width:780px; margin:0 auto; padding-left:max(16px, env(safe-area-inset-left,0px)); padding-right:max(16px, env(safe-area-inset-right,0px)); padding-block:40px 80px; }
h1,h2,h3 { font-family:var(--display); font-weight:600; text-wrap:balance; line-height:1.2; }
h1 { font-size:clamp(1.9rem,4.6vw,2.6rem); margin:0 0 16px; }
h2 { font-size:1.55rem; margin:44px 0 12px; padding-top:12px; border-top:1px solid var(--line); }
h3 { font-size:1.2rem; margin:28px 0 8px; }
p, li { max-width:70ch; }
p { margin:0 0 14px; }
ul, ol { padding-left:22px; margin:0 0 14px; }
li { margin:4px 0; }
a { color:var(--accent); }
hr { border:0; border-top:1px solid var(--line); margin:28px 0; }
em { font-style:italic; }
code { font-family:var(--mono); font-size:0.88em; }
.tw { overflow-x:auto; margin:8px 0 18px; }
table { border-collapse:collapse; font-size:0.92rem; min-width:100%; }
th, td { text-align:left; vertical-align:top; padding:7px 10px; border-bottom:1px solid var(--line); }
th { background:var(--tablehead); font-weight:600; }
figure { margin:22px 0 26px; }
figure img { display:block; width:100%; height:auto; background:var(--plate); border:1px solid var(--line); border-radius:8px; padding:8px; box-sizing:border-box; }
figcaption { font-size:0.9rem; color:var(--muted); margin-top:8px; max-width:70ch; }
main > p:first-of-type em { color:var(--muted); }
</style>
<main>
"""
out = Path(sys.argv[1]); out.write_text(head + body + "\n</main>\n")
print(out, out.stat().st_size // 1024, "KB", body.count("<figure>"), "figures")

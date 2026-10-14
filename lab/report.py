"""Render out/citations.json from citation_sweep.py.

ASCII-safe: Windows consoles default to cp1252 and OpenAlex titles are full of
en-dashes and smart quotes, which raise UnicodeEncodeError on print.

Usage:
    python report.py            # works citing 2+ anchors (default)
    python report.py recent     # 2025 onward
    python report.py all
    python report.py shear      # any keyword match against the title
"""
import json
import os
import re
import sys

P = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "citations.json")

if not os.path.exists(P):
    sys.exit(f"no {P} -- run citation_sweep.py first")

rows = json.load(open(P, encoding="utf-8"))

_SUBS = {"‐": "-", "‑": "-", "‒": "-", "–": "-",
         "—": "-", "‘": "'", "’": "'", "“": '"',
         "”": '"', " ": " ", "−": "-"}


def safe(s, n=None):
    s = s or ""
    for a, b in _SUBS.items():
        s = s.replace(a, b)
    s = re.sub(r"\s+", " ", s.encode("ascii", "replace").decode("ascii")).strip()
    return s[:n] if n else s


mode = sys.argv[1] if len(sys.argv) > 1 else "multi"
if mode == "multi":
    sel, title = [r for r in rows if r["n_anchors"] >= 2], "MULTI-ANCHOR (2+)"
elif mode == "all":
    sel, title = rows, "ALL"
elif mode == "recent":
    sel, title = [r for r in rows if r["date"] >= "2025-01-01"], "2025 ONWARD"
else:
    sel = [r for r in rows if mode.lower() in (r["title"] or "").lower()]
    title = f"KEYWORD: {mode}"

sel.sort(key=lambda r: (r["date"] or str(r["year"]), r["n_anchors"]), reverse=True)

print("=" * 78)
print(f"{title} -- {len(sel)} of {len(rows)} works, newest first")
print("=" * 78)
for r in sel:
    print(f"[{r['date'] or r['year']}] ({r['n_anchors']}x: {','.join(r['cites'])}) "
          f"cb={r['cited_by']}")
    print(f"   {safe(r['title'], 106)}")
    print(f"   {safe(r['venue'], 56)} | {r['doi']}")

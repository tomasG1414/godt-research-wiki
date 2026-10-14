#!/usr/bin/env python3
"""Mechanical wiki lint. Stdlib only. Exit 0 clean, 1 if issues.

Checks:
  - broken [[wikilinks]] (target page does not exist)
  - orphan pages (no inbound links from anywhere, index included)
  - duplicate page slugs (two files with the same stem in different folders)
  - index drift (pages missing from index.md)

Judgment checks (contradictions, staleness, missing concepts) are the agent's
job — see CLAUDE.md, Lint.
"""
import re
import sys
from collections import defaultdict
from pathlib import Path

# [[target]] / [[target|alias]] / [[target#heading]] -> capture target only
WIKILINK = re.compile(r"\[\[([^\]\|#]+)")
SPECIAL = {"index"}


def main():
    root = Path(__file__).resolve().parent.parent / "wiki"
    files = sorted(root.rglob("*.md"))
    if not files:
        sys.exit(f"no .md files under {root}")

    by_stem = defaultdict(list)
    for f in files:
        by_stem[f.stem].append(f)
    pages = set(by_stem)

    duplicates = {s: fs for s, fs in by_stem.items() if len(fs) > 1}

    inbound = {s: 0 for s in pages}
    broken = []
    indexed = set()
    for f in files:
        text = f.read_text(errors="ignore")
        for m in WIKILINK.finditer(text):
            target = Path(m.group(1).strip()).stem  # tolerate path-style links
            if target in pages:
                if target != f.stem:
                    inbound[target] += 1
                if f.stem == "index":
                    indexed.add(target)
            else:
                broken.append((f.relative_to(root), m.group(1).strip()))

    orphans = sorted(s for s in pages if s not in SPECIAL and inbound[s] == 0)
    unindexed = sorted(s for s in pages if s not in SPECIAL and s not in indexed)

    issues = 0

    def report(title, items, fmt=str):
        nonlocal issues
        if items:
            issues += len(items)
            print(f"{title} ({len(items)}):")
            for it in items:
                print(f"  {fmt(it)}")
            print()

    report("Broken wikilinks", broken, lambda x: f"{x[0]} -> [[{x[1]}]]")
    report("Orphan pages (no inbound links)", orphans)
    report("Duplicate slugs", sorted(duplicates.items()),
           lambda x: f"{x[0]}: " + ", ".join(str(p.relative_to(root)) for p in x[1]))
    report("Pages missing from index.md", unindexed)

    if issues == 0:
        print("clean")
    else:
        print(f"{issues} issue(s)")
    sys.exit(0 if issues == 0 else 1)


if __name__ == "__main__":
    main()

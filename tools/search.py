#!/usr/bin/env python3
"""BM25 search over the wiki. Stdlib only.

Usage:
    python3 tools/search.py "ewald cap rytov" [-k 8] [--dir PATH]
"""
import argparse
import math
import re
import sys
from collections import Counter
from pathlib import Path

TOKEN = re.compile(r"[a-z0-9]+")
K1, B = 1.5, 0.75


def tokenize(text: str):
    return TOKEN.findall(text.lower())


def main():
    ap = argparse.ArgumentParser(description="BM25 search over wiki markdown pages")
    ap.add_argument("query", nargs="+", help="query terms")
    ap.add_argument("-k", type=int, default=8, help="max results (default 8)")
    ap.add_argument("--dir", default=None, help="folder to search (default: wiki/)")
    args = ap.parse_args()

    root = Path(args.dir) if args.dir else Path(__file__).resolve().parent.parent / "wiki"
    files = sorted(root.rglob("*.md"))
    if not files:
        sys.exit(f"no .md files under {root}")

    texts = [f.read_text(errors="ignore") for f in files]
    docs = [Counter(tokenize(t)) for t in texts]
    n_docs = len(docs)
    avgdl = sum(sum(d.values()) for d in docs) / n_docs

    query = tokenize(" ".join(args.query))
    if not query:
        sys.exit("empty query")
    df = {t: sum(1 for d in docs if t in d) for t in set(query)}

    scored = []
    for i, d in enumerate(docs):
        dl = sum(d.values()) or 1
        s = 0.0
        for t in query:
            if t not in d or df[t] == 0:
                continue
            idf = math.log(1 + (n_docs - df[t] + 0.5) / (df[t] + 0.5))
            tf = d[t]
            s += idf * tf * (K1 + 1) / (tf + K1 * (1 - B + B * dl / avgdl))
        if s > 0:
            scored.append((s, i))
    scored.sort(reverse=True)

    if not scored:
        print("no matches")
        return

    base = root.parent
    for s, i in scored[: args.k]:
        try:
            rel = files[i].relative_to(base)
        except ValueError:
            rel = files[i]
        print(f"{s:6.2f}  {rel}")
        shown = 0
        for line in texts[i].splitlines():
            stripped = line.strip()
            if stripped and any(t in stripped.lower() for t in query):
                print(f"        {stripped[:120]}")
                shown += 1
                if shown == 2:
                    break


if __name__ == "__main__":
    main()

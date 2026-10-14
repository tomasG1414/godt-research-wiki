"""Forward citation sweep over OpenAlex for the project's load-bearing anchors.

Resolves each anchor to an OpenAlex work ID, pulls everything citing it since
FROM_YEAR, and aggregates so works citing MULTIPLE anchors surface first --
those sit in the same intellectual neighbourhood as the project.

Backward citations tell you where a paper came from; this tells you who picked
it up. That is where recent work lives.

Writes: out/citations.json   (render it with report.py)
Wiki:   wiki/notes/forward-citation-sweep-2026-07.md

Usage:
    python citation_sweep.py
    OPENALEX_MAILTO=you@example.com python citation_sweep.py

OpenAlex asks for a contact address to put you in the polite pool. Set
OPENALEX_MAILTO or edit the default below.
"""
import json
import os
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict

MAILTO = os.environ.get("OPENALEX_MAILTO", "research-wiki@example.com")
BASE = "https://api.openalex.org"
FROM_YEAR = 2024

# (key, DOI or None, fallback title search)
# NOTE: chen1998 has never resolved -- its DOI lookup returns nothing and the
# title fallback matches an unrelated paper, so it is disabled rather than
# silently polluting results. Its forward citations remain unswept.
ANCHORS = [
    ("wolf1969", "10.1016/0030-4018(69)90052-2",
     "Three-dimensional structure determination of semi-transparent objects "
     "from holographic data"),
    ("brenner1993", "10.1364/AO.32.004984",
     "Light propagation through microlenses: a new simulation method"),
    ("vishnyakov2016", "10.1364/OL.41.003037",
     "Differential interference contrast tomography"),
    ("krauze2016", "10.1364/OE.24.004924",
     "Generalized total variation iterative constraint strategy in limited "
     "angle optical diffraction tomography"),
    ("nguyen2017", "10.1038/s41467-017-00190-7",
     "Gradient light interference microscopy for 3D imaging of unlabeled specimens"),
    ("ziemczonok2019", "10.1038/s41598-019-55330-4",
     "3D-printed biological cell phantom for testing 3D quantitative phase "
     "imaging systems"),
    ("godt2025", "10.1038/s42005-025-02416-3",
     "Gradient optical diffraction tomography"),
    ("osiflex2025", "10.1088/2515-7647/ae3504",
     "OSI-flex optimization-based shearing interferometry joint phase shear"),
    # ("chen1998", "10.1364/AO.37.002996", "Validity of diffraction tomography "
    #  "based on the first Born and the first Rytov approximations"),  # unresolved
]


def out_path(fname):
    d = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, fname)


def get(url, tries=3):
    req = urllib.request.Request(
        url, headers={"User-Agent": f"research-wiki-sweep (mailto:{MAILTO})"})
    for a in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.load(r)
        except Exception as e:
            if a == tries - 1:
                print(f"    FAILED {url}\n      {e}", file=sys.stderr)
                return None
            time.sleep(2 * (a + 1))
    return None


def resolve(doi, title):
    if doi:
        d = get(f"{BASE}/works/doi:{urllib.parse.quote(doi)}?mailto={MAILTO}")
        if d and d.get("id"):
            return (d["id"].rsplit("/", 1)[-1], d.get("display_name"),
                    d.get("publication_year"), d.get("cited_by_count", 0))
    d = get(f"{BASE}/works?search={urllib.parse.quote(title)}"
            f"&per-page=5&mailto={MAILTO}")
    if d and d.get("results"):
        t = d["results"][0]
        return (t["id"].rsplit("/", 1)[-1], t.get("display_name"),
                t.get("publication_year"), t.get("cited_by_count", 0))
    return None


def citing(work_id, from_year):
    out, cursor = [], "*"
    while cursor:
        d = get(f"{BASE}/works?filter=cites:{work_id},"
                f"from_publication_date:{from_year}-01-01"
                f"&per-page=200&cursor={urllib.parse.quote(cursor)}&mailto={MAILTO}")
        if not d:
            break
        out.extend(d.get("results", []))
        cursor = d.get("meta", {}).get("next_cursor")
        time.sleep(0.15)
    return out


def main():
    print("=" * 78 + "\nRESOLVING ANCHORS\n" + "=" * 78)
    resolved = {}
    for key, doi, title in ANCHORS:
        r = resolve(doi, title)
        if not r:
            print(f"  {key:15s} UNRESOLVED")
            continue
        wid, name, yr, tot = r
        resolved[key] = wid
        nm = (name or "").encode("ascii", "replace").decode("ascii")
        print(f"  {key:15s} {wid:12s} {yr}  cites={tot:<6d} {nm[:46]}")
        time.sleep(0.15)

    print("\n" + "=" * 78 + f"\nFORWARD CITATIONS ({FROM_YEAR}+)\n" + "=" * 78)
    universe, hits = {}, defaultdict(set)
    for key, wid in resolved.items():
        works = citing(wid, FROM_YEAR)
        print(f"  {key:15s} {len(works):4d} citing works")
        for w in works:
            cid = w["id"].rsplit("/", 1)[-1]
            universe[cid] = w
            hits[cid].add(key)

    rows = []
    for cid, w in universe.items():
        src = ((w.get("primary_location") or {}).get("source")
               or {}).get("display_name", "")
        rows.append({
            "id": cid, "title": w.get("display_name") or "",
            "year": w.get("publication_year") or 0,
            "date": w.get("publication_date") or "",
            "venue": src or "",
            "doi": (w.get("doi") or "").replace("https://doi.org/", ""),
            "cites": sorted(hits[cid]), "n_anchors": len(hits[cid]),
            "cited_by": w.get("cited_by_count", 0), "type": w.get("type") or "",
        })
    rows.sort(key=lambda r: (r["n_anchors"], r["date"]), reverse=True)

    p = out_path("citations.json")
    with open(p, "w", encoding="utf-8") as f:
        json.dump(rows, f, indent=1)
    print(f"\n  {len(rows)} unique citing works; "
          f"{sum(1 for r in rows if r['n_anchors'] >= 2)} cite 2+ anchors")
    print(f"  wrote {p}")
    print("  render with:  python report.py multi|recent|all|<keyword>")


if __name__ == "__main__":
    main()

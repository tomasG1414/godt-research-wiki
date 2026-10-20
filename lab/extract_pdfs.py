"""Extract text from the PDFs in raw/ so they can be read during an ingest.

Uses PyMuPDF for text-layer extraction, with no poppler dependency.
Trade-off worth remembering -- this reads TEXT ONLY, so figures and anything
living only in an image are invisible. Figure captions usually survive.

Writes: out/extracted/<name>.txt, page-delimited.

Usage:
    python extract_pdfs.py                    # default: raw/godt
    python extract_pdfs.py raw/other/folder
"""
import os
import sys

try:
    import fitz  # PyMuPDF
except ImportError:
    sys.exit("PyMuPDF missing.  pip install -r lab/requirements-ingest.txt")

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
DEFAULT_SRC = os.path.join(REPO, "raw", "godt")


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_SRC
    if not os.path.isabs(src):
        src = os.path.join(REPO, src)
    if not os.path.isdir(src):
        sys.exit(f"not a directory: {src}")

    out = os.path.join(HERE, "out", "extracted")
    os.makedirs(out, exist_ok=True)

    pdfs = sorted(f for f in os.listdir(src) if f.lower().endswith(".pdf"))
    if not pdfs:
        sys.exit(f"no PDFs in {src}")
    print(f"{len(pdfs)} PDFs in {src}")

    ok = 0
    for fname in pdfs:
        try:
            doc = fitz.open(os.path.join(src, fname))
            parts = []
            for i, page in enumerate(doc):
                parts.append(f"\n\n===== PAGE {i+1}/{len(doc)} =====\n\n")
                parts.append(page.get_text())
            text = "".join(parts)
            dest = os.path.join(out, fname.rsplit(".", 1)[0] + ".txt")
            with open(dest, "w", encoding="utf-8") as f:
                f.write(text)
            print(f"  OK   {fname:46s} {len(doc):3d}p {len(text):7d} chars")
            ok += 1
        except Exception as e:
            print(f"  FAIL {fname}: {e}")
    print(f"\n{ok}/{len(pdfs)} extracted -> {out}")


if __name__ == "__main__":
    main()

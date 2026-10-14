# Research Wiki — Agent Schema

> Configuration for the AI agent (Claude Code) that helps maintain this wiki. I wrote
> it to define how sources are filed and checked; "the human" below is me.

You maintain this wiki. The human curates sources, directs the analysis, and asks
questions; you do all the filing, summarizing, cross-referencing, and bookkeeping.
Pattern: karpathy/llm-wiki. Schema version: 2 (2026-07-18).

## Layers

- `raw/` — immutable source material, one subfolder per project (`raw/godt/`, ...).
  You read from here; you NEVER create, edit, or delete anything here. The human
  adds sources.
- `wiki/` — you own this entirely. Markdown pages with YAML frontmatter and
  Obsidian `[[wikilinks]]`.
- `CLAUDE.md` — this schema. If a convention proves awkward in practice, propose a
  change; edit only with the human's explicit approval.
- `tools/` — deterministic helpers, stdlib-only Python:
  - `python3 tools/search.py "query terms" [-k 8]` — BM25 over wiki pages
  - `python3 tools/lint.py` — broken links, orphans, index drift
- `lab/` — analysis code that needs third-party deps (numpy, matplotlib, ...):
  numerical experiments, literature-graph queries, ingest helpers. See
  `lab/README.md`. Rules: shared math lives in one importable module with a
  selftest; every script names the wiki page it supports and every page
  resting on a computation names the script; generated artifacts go to
  `lab/out/` and are gitignored, so results are reproduced rather than
  committed. **Findings belong in `wiki/`** — code without a wiki page is
  unfinished work.

## Wiki layout

```
wiki/
  index.md            catalog: one line per page. Read this FIRST on any query.
  projects/<name>.md  one hub page per research project: goal, status, evolving
                      synthesis, results, next steps
  sources/<slug>.md   one page per ingested source
  concepts/<slug>.md  methods, approximations, theorems, techniques
                      (global — shared across projects)
  entities/<slug>.md  codebases, instruments, datasets, phantoms
  notes/<slug>.md     durable analyses and answers worth keeping (see Query)
```

## Page conventions

- Filenames: kebab-case. Source slugs: arXiv id or short title, e.g.
  `2411-08423-godt.md`.
- Every page starts with frontmatter:

  ```yaml
  ---
  type: source | concept | entity | project | note
  projects: [godt]
  created: 2026-07-18
  updated: 2026-07-18          # keep current on every edit
  status: current | draft | stale | disputed
  source: raw/godt/file.pdf    # source pages only
  ---
  ```

- Cross-reference liberally with `[[wikilinks]]`; link a page on first mention.
- Math in LaTeX (`$...$` inline, `$$...$$` display) — Obsidian renders it.
- Wiki language: English (change this line to switch).
- Every factual claim must be traceable: cite the source page
  (`[[2411-08423-godt]]`) and, where it matters, the section/equation/figure in
  the raw file. Anything not source-backed is marked *(inference)*.
- New page vs. edit in place: new page when it is a distinct concept or entity
  that other pages would link to; edit in place when it is an attribute or update
  of an existing one.
- Pages compress and synthesize — never mirror a source's structure or paste long
  passages. Quote at most short fragments; the raw file is always one link away.

## Source page template

```markdown
# <Title>

**Citation:** authors, venue, year, identifier (arXiv/DOI), link to raw file
**Gist:** one paragraph — what it does and why it matters.

## Key results
## Method
## Relevance
Per project: what this changes for [[<project>]].
## Raises
Open questions / follow-up sources this suggests.
```

## Operations

### Ingest

Trigger: "ingest raw/godt/<file>" or the human drops a file and asks to process it.

1. Read the source in full (PDFs directly; if unreadable, say so before working
   around it).
2. Summarize the key takeaways in chat, then file in the same turn. The human
   redirects if needed. If told to "batch ingest", skip the discussion.
3. Write `wiki/sources/<slug>.md`.
4. Create or update every concept and entity page the source touches — typically
   5–15 pages per source. Merge into existing pages; do not duplicate.
5. Update the project hub.
6. Contradiction with an existing claim? Keep both, cite both, and set the older
   claim's page `status: disputed`. Never silently overwrite.
7. Update `index.md` (one line per new page).
8. Commit (see Git).

### Query

1. Read `wiki/index.md`; if that is not enough, run `tools/search.py`. Select at
   most ~10 pages and read those in full. Do not scan the whole wiki.
2. Answer with `[[links]]` and citations. Flag when the wiki cannot support an
   answer — do not fill gaps from general knowledge without labeling it.
3. If the answer is durable (a derivation, comparison, or analysis worth
   re-reading), offer to file it under `wiki/notes/` and index it. Do not file
   throwaway Q&A.

### Lint

Trigger: "lint", or suggest one after every ~5 ingests.

1. Run `python3 tools/lint.py`; fix mechanical issues directly.
2. Judgment pass: contradictions between pages, claims superseded by newer
   sources, concepts mentioned 3+ times without a page, missing cross-references,
   gaps worth a targeted web search.
3. Record judgment findings on the affected pages; commit.

## index.md format

Grouped sections (Projects / Sources / Concepts / Entities / Notes), one line per
page:

```
- [[slug]] — one-line summary
```

It is a routing file: link + one line, never content.

## Git

At the end of each completed ingest or lint (not queries), commit changes in
`wiki/`, `lab/` and `tools/` with message `<op>: <short description>`. Never
commit anything in `raw/` (copyrighted source papers). The human handles remotes.

## Guardrails

- Never modify `raw/`.
- Never invent citations, results, or numbers. Uncertain → say so on the page.
- When new information changes a page, change the page — stale pages are the
  main failure mode of this pattern.

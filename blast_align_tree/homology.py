"""Homology report: how weak the weakest retained BLAST hits are.

The number of hits kept per search is fixed by -n (-max_target_seqs) rather
than by an e-value cutoff, because a sensible cutoff depends on the gene family,
sequence length and database size. Instead of choosing one, every run reports,
for each (query, database) search, the worst hit it kept: its e-value, bit
score, identity, similarity and query coverage. Whether -n truncated the search
is stated alongside, so a reader can judge whether weaker homologs may have been
left out.
"""

from __future__ import annotations

import csv
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

# Columns BLAST writes per hit (-outfmt "6 ..."), in this order.
BLAST_HIT_FIELDS = ["sseqid", "pident", "ppos", "length", "qcovhsp", "evalue", "bitscore"]
BLAST_OUTFMT = "6 " + " ".join(BLAST_HIT_FIELDS)
HITS_SUFFIX = ".hits.tsv"

# BLAST+'s default -evalue; no other e-value threshold is applied.
BLAST_DEFAULT_EVALUE = 10

LOG_COLUMNS = [
    "query", "database", "max_target_seqs", "hits_returned", "limited_by_n",
    "worst_hit", "worst_hit_identifier", "worst_evalue", "worst_bitscore",
    "worst_pct_identity", "worst_pct_similarity", "worst_aln_length",
    "worst_query_coverage", "best_evalue", "best_pct_identity",
]

COLUMN_DESCRIPTION = {
    "limited_by_n": "yes if hits_returned reached max_target_seqs, so weaker "
                    "homologs may exist beyond the -n cap; no if every hit "
                    f"with e-value <= {BLAST_DEFAULT_EVALUE} was kept",
    "worst_hit": "the retained hit with the highest e-value (ties: lowest bit "
                 "score), as named in the BLAST database",
    "worst_hit_identifier": "its label in the tree; '-' if dropped by de-duplication",
    "worst_pct_identity": "BLAST pident: identical positions over the HSP",
    "worst_pct_similarity": "BLAST ppos: positive-scoring positions over the HSP",
    "worst_query_coverage": "BLAST qcovhsp: percent of the query covered by the HSP",
}


def read_hits(path: Path) -> List[Dict[str, object]]:
    """Parse one search's per-hit table written with BLAST_OUTFMT."""
    hits: List[Dict[str, object]] = []
    if not path.exists():
        return hits
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) != len(BLAST_HIT_FIELDS):
            continue
        hit: Dict[str, object] = dict(zip(BLAST_HIT_FIELDS, parts))
        for key in ("pident", "ppos", "evalue", "bitscore"):
            hit[key] = float(parts[BLAST_HIT_FIELDS.index(key)])
        for key in ("length", "qcovhsp"):
            hit[key] = int(float(parts[BLAST_HIT_FIELDS.index(key)]))
        hits.append(hit)
    return hits


def _worst(hits: Sequence[Dict[str, object]]) -> Dict[str, object]:
    return max(hits, key=lambda h: (h["evalue"], -h["bitscore"]))


def _best(hits: Sequence[Dict[str, object]]) -> Dict[str, object]:
    return min(hits, key=lambda h: (h["evalue"], -h["bitscore"]))


def _fmt_evalue(x: float) -> str:
    return "0" if x == 0 else f"{x:.2e}"


def build_row(query: str, database: str, max_targets: str,
              hits: Sequence[Dict[str, object]],
              final_ids: Dict[Tuple[str, str], str]) -> Dict[str, object]:
    """One report row for a (query, database) search."""
    row: Dict[str, object] = {
        "query": query, "database": database, "max_target_seqs": max_targets,
        "hits_returned": len(hits),
        "limited_by_n": "yes" if len(hits) >= int(max_targets) else "no",
    }
    if not hits:
        for col in LOG_COLUMNS[5:]:
            row[col] = "-"
        return row
    worst, best = _worst(hits), _best(hits)
    row.update({
        "worst_hit": worst["sseqid"],
        "worst_hit_identifier": final_ids.get((database, str(worst["sseqid"])), "-"),
        "worst_evalue": _fmt_evalue(worst["evalue"]),
        "worst_bitscore": f"{worst['bitscore']:g}",
        "worst_pct_identity": f"{worst['pident']:.1f}",
        "worst_pct_similarity": f"{worst['ppos']:.1f}",
        "worst_aln_length": worst["length"],
        "worst_query_coverage": worst["qcovhsp"],
        "best_evalue": _fmt_evalue(best["evalue"]),
        "best_pct_identity": f"{best['pident']:.1f}",
    })
    return row


def write_log(path: Path, rows: Sequence[Dict[str, object]], *, entry: str,
              blast_type: str) -> Path:
    """Write the per-search homology report with a self-describing header."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        fh.write("# blast-align-tree homology report\n")
        fh.write(f"# entry: {entry}\n")
        fh.write(f"# generated: {datetime.now().isoformat(timespec='seconds')}\n")
        fh.write(f"# blast_type: {blast_type}\n")
        fh.write("# hits per search are capped by -n (-max_target_seqs, one HSP per "
                 f"subject); no e-value cutoff beyond BLAST's default ({BLAST_DEFAULT_EVALUE})\n")
        fh.write("# one row per (query, database) search; worst_* describe the weakest "
                 "hit retained, best_* the strongest\n")
        for col, desc in COLUMN_DESCRIPTION.items():
            fh.write(f"# {col}: {desc}\n")
        writer = csv.DictWriter(fh, fieldnames=LOG_COLUMNS, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


def print_report(rows: Sequence[Dict[str, object]], log_path: Optional[Path]) -> None:
    """End-of-run summary: the weakest hit kept, and which searches hit the cap."""
    with_hits = [r for r in rows if r["hits_returned"]]
    print(f"\n  [homology] {len(rows)} BLAST searches; worst retained hit per search "
          f"in {log_path.name if log_path else 'the homology report'}")
    if with_hits:
        weakest = max(with_hits, key=lambda r: float(r["worst_evalue"]))
        print(f"    Weakest hit kept: {weakest['worst_hit']} ({weakest['query']} vs "
              f"{weakest['database']}): e-value {weakest['worst_evalue']}, "
              f"{weakest['worst_pct_identity']}% identity, "
              f"{weakest['worst_pct_similarity']}% similarity")
    capped = [r for r in rows if r["limited_by_n"] == "yes"]
    if capped:
        print(f"    {len(capped)} of {len(rows)} searches reached the -n cap; "
              "weaker homologs may exist beyond it")
    if log_path:
        print(f"    Report: {log_path}")

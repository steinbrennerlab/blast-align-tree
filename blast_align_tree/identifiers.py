#!/usr/bin/env python3
"""
Deterministic, auditable handling of identifier collisions after header parsing.

Header parsing (``-hdr gene:``, ``-hdr locus=``, ...) turns a source record's
description into a short identifier used as the tree tip label. That mapping is
not injective: two isoforms of one locus share a gene symbol, and two genomes
can independently use the same symbol. Collapsing those to "first one wins"
throws away records with no record of what was lost.

This module makes every collision explicit:

  * **isoform / duplicated locus** - one database, one identifier, several
    distinct source records. The longest amino-acid sequence is retained; ties
    break on the lexicographically smallest source ID. The loser is dropped and
    logged with the ID of the record that replaced it.
  * **cross-query** - the same shape as above, except no single query hit more
    than one of the records, so it only becomes visible when queries are merged.
    Same rule, same log, no separate prompt.
  * **cross-database** - one identifier appearing in more than one genome.
    Dropping here would delete a species rather than a redundant isoform, so
    every member survives with a genome tag appended.

Two identical identifiers backed by the *same* source record are not a
collision; that is ordinary overlap between queries and is reported separately.

Everything here is pure functions over paths and dataclasses, so the policy can
be tested without running BLAST.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

try:
    from Bio import SeqIO
except ImportError:  # pragma: no cover - cli.py reports the missing dependency
    SeqIO = None


# Collision kinds
ISOFORM = "isoform"
CROSS_QUERY = "cross_query"
CROSS_DATABASE = "cross_database"

# Log stages
STAGE_BY_KIND = {
    ISOFORM: "within_query",
    CROSS_QUERY: "cross_query",
    CROSS_DATABASE: "cross_database",
}


# -----------------------
# Header token parsing
# -----------------------

def parse_header_token(description: str, headerword: str, fallback_id: str, suffix: str = "") -> str:
    """
    If headerword == 'id' -> return fallback_id (record.id).
    Else: find substring after the first occurrence of headerword and take until next space.
    If headerword not found, return fallback_id.
    If suffix is provided, strip it from the end of the token.
    """
    if headerword == "id":
        token = fallback_id
    elif headerword not in description:
        token = fallback_id
    else:
        # Split on the headerword and take the part immediately following, up to first space
        try:
            part = description.split(headerword, 1)[1]
            token = part.split(" ", 1)[0]
        except Exception:
            token = fallback_id
    if suffix and token.endswith(suffix):
        token = token[:-len(suffix)]
    return token


# -----------------------
# Data model
# -----------------------

@dataclass(frozen=True)
class RecordRef:
    """One source record as hit by one query in one database."""
    query: str
    db: str
    source_id: str
    description: str
    aa_len: int
    token: str


@dataclass
class Collision:
    """A single identifier claimed by more than one distinct source record."""
    kind: str
    token: str
    db: Optional[str]              # None for CROSS_DATABASE (spans databases)
    members: List[RecordRef]       # one per distinct source record, best first
    queries: Tuple[str, ...]       # queries implicated in this collision

    @property
    def winner(self) -> RecordRef:
        return self.members[0]

    @property
    def losers(self) -> List[RecordRef]:
        return self.members[1:]


@dataclass
class Resolution:
    """The final identifier assigned to every record, and what was dropped."""
    final_ids: Dict[Tuple[str, str], str] = field(default_factory=dict)  # (db, source_id) -> identifier
    dropped: Set[Tuple[str, str]] = field(default_factory=set)           # (db, source_id)
    collisions: List[Collision] = field(default_factory=list)
    declined_tokens: Set[Tuple[str, str]] = field(default_factory=set)   # (db, token) kept via suffixes

    def submap(self, db: str) -> Dict[str, str]:
        """source_id -> final identifier, for one database."""
        return {sid: fid for (d, sid), fid in self.final_ids.items() if d == db}

    def dropped_in(self, db: str) -> Set[str]:
        return {sid for (d, sid) in self.dropped if d == db}


# -----------------------
# Index construction
# -----------------------

def aa_source_path(entry_dir: Path, query: str, db_label: str, bt: str, blast_type: str) -> Path:
    """The pre-parse FASTA carrying amino-acid sequences for one (query, db)."""
    if blast_type == "tblastn":
        return entry_dir / f"{query}.{db_label}.seq.{bt}.blastdb.translate.fa"
    return entry_dir / f"{query}.{db_label}.seq.{bt}.blastdb.fa"


def build_index(
    entry_dir: Path,
    queries: Sequence[str],
    databases: Sequence[str],
    hdr_rules_by_db: Dict[str, Tuple[str, str]],
    blast_type: str,
    bt: str,
    db_label_fn,
) -> List[RecordRef]:
    """
    Read every (query, db) pre-parse FASTA and recompute the parsed identifier
    alongside the source ID it came from.

    The parsed FASTAs written during the BLAST loop keep only the token, so the
    source ID and original description have to be recovered here - that pairing
    is the whole basis for telling a real collision from benign overlap.
    """
    index: List[RecordRef] = []
    for db in databases:
        dbl = db_label_fn(db)
        headerword, suffix = hdr_rules_by_db[db]
        for query in queries:
            fp = aa_source_path(entry_dir, query, dbl, bt, blast_type)
            if not fp.is_file():
                continue
            for rec in SeqIO.parse(str(fp), "fasta"):
                token = parse_header_token(rec.description, headerword, rec.id, suffix)
                index.append(RecordRef(
                    query=query,
                    db=db,
                    source_id=rec.id,
                    description=rec.description,
                    aa_len=len(rec.seq),
                    token=token,
                ))
    return index


# -----------------------
# Detection
# -----------------------

def _rank_key(ref: RecordRef) -> Tuple[int, str]:
    """Longest amino-acid sequence first; ties break on smallest source ID."""
    return (-ref.aa_len, ref.source_id)


def _distinct_by_source(refs: Iterable[RecordRef]) -> List[RecordRef]:
    """One RecordRef per source ID, ranked best first."""
    best: Dict[str, RecordRef] = {}
    for ref in refs:
        best.setdefault(ref.source_id, ref)
    return sorted(best.values(), key=_rank_key)


def detect(index: Sequence[RecordRef]) -> List[Collision]:
    """
    Find every identifier backed by more than one distinct source record.

    Runs before any user confirmation, so it must not depend on the outcome.
    """
    by_db_token: Dict[Tuple[str, str], List[RecordRef]] = defaultdict(list)
    for ref in index:
        by_db_token[(ref.db, ref.token)].append(ref)

    collisions: List[Collision] = []

    # Within a database: isoforms, duplicated loci.
    for (db, token), refs in sorted(by_db_token.items()):
        members = _distinct_by_source(refs)
        if len(members) < 2:
            continue

        # A query that hit two or more of these records collides with itself and
        # is worth confirming. If every query hit only one, the clash only
        # appears once queries are merged.
        per_query: Dict[str, Set[str]] = defaultdict(set)
        for ref in refs:
            per_query[ref.query].add(ref.source_id)
        self_colliding = sorted(q for q, ids in per_query.items() if len(ids) > 1)

        if self_colliding:
            kind, queries = ISOFORM, tuple(self_colliding)
        else:
            kind, queries = CROSS_QUERY, tuple(sorted(per_query))

        collisions.append(Collision(kind=kind, token=token, db=db, members=members, queries=queries))

    # Across databases: the same symbol used by different genomes.
    dbs_by_token: Dict[str, Set[str]] = defaultdict(set)
    for (db, token) in by_db_token:
        dbs_by_token[token].add(db)

    for token, dbs in sorted(dbs_by_token.items()):
        if len(dbs) < 2:
            continue
        members = [ref for db in sorted(dbs) for ref in _distinct_by_source(by_db_token[(db, token)])]
        queries = tuple(sorted({r.query for db in dbs for r in by_db_token[(db, token)]}))
        collisions.append(Collision(kind=CROSS_DATABASE, token=token, db=None,
                                    members=members, queries=queries))

    return collisions


def collisions_by_query(collisions: Sequence[Collision]) -> Dict[str, List[Collision]]:
    """Within-query collisions grouped by the query that needs to confirm them."""
    grouped: Dict[str, List[Collision]] = defaultdict(list)
    for col in collisions:
        if col.kind != ISOFORM:
            continue
        for query in col.queries:
            grouped[query].append(col)
    return grouped


# -----------------------
# Genome tags
# -----------------------

def genome_tag(db: str) -> str:
    """Short, readable label for a database: 'NbLab360.v103.gff3.CDS.fasta' -> 'NbLab360'."""
    stem = Path(db).name
    for ext in (".fasta", ".fa", ".fna", ".faa", ".pep"):
        if stem.lower().endswith(ext):
            stem = stem[: -len(ext)]
            break
    return stem.split(".", 1)[0] or stem


def _tags_for(dbs: Sequence[str]) -> Dict[str, str]:
    """
    Genome tags guaranteed distinct within one collision. Falls back to the full
    filename stem, then to the flattened path, rather than emitting two records
    that would collide all over again.
    """
    tags = {db: genome_tag(db) for db in dbs}
    if len(set(tags.values())) == len(dbs):
        return tags

    tags = {db: Path(db).name.rsplit(".", 1)[0] for db in dbs}
    if len(set(tags.values())) == len(dbs):
        return tags

    return {db: db.replace("/", ".").replace("\\", ".") for db in dbs}


# -----------------------
# Resolution
# -----------------------

def resolve(
    index: Sequence[RecordRef],
    collisions: Sequence[Collision],
    declined: Iterable[Tuple[str, str]] = (),
) -> Resolution:
    """
    Assign a final identifier to every record.

    ``declined`` holds (db, token) pairs the user refused to deduplicate; those
    keep every member, disambiguated by rank suffix. Resolution is computed
    globally per (db, token) so that every query agrees on the same outcome.
    """
    declined_set = set(declined)
    res = Resolution(collisions=list(collisions), declined_tokens=declined_set)

    # Baseline: every record keeps its parsed token.
    for ref in index:
        res.final_ids[(ref.db, ref.source_id)] = ref.token

    # Stage 1 - within a database, pick a winner (or suffix everything).
    suffixed: Dict[Tuple[str, str], int] = {}
    for col in collisions:
        if col.kind == CROSS_DATABASE:
            continue
        key = (col.db, col.token)
        if key in declined_set:
            for rank, ref in enumerate(col.members, start=1):
                suffixed[(ref.db, ref.source_id)] = rank
        else:
            for ref in col.losers:
                res.dropped.add((ref.db, ref.source_id))
                res.final_ids.pop((ref.db, ref.source_id), None)

    # Stage 2 - across databases, tag survivors with their genome.
    for col in collisions:
        if col.kind != CROSS_DATABASE:
            continue
        dbs = sorted({ref.db for ref in col.members})
        tags = _tags_for(dbs)
        for (db, source_id), token in list(res.final_ids.items()):
            if token == col.token and db in tags:
                res.final_ids[(db, source_id)] = f"{col.token}_{tags[db]}"

    # Rank suffixes apply last so a declined cross-database token reads
    # <token>_<genome>__<rank> rather than interleaving the two markers.
    for key, rank in suffixed.items():
        if key in res.final_ids:
            res.final_ids[key] = f"{res.final_ids[key]}__{rank}"

    return res


# -----------------------
# Overlap between queries
# -----------------------

@dataclass(frozen=True)
class QueryOverlap:
    db: str
    query_a: str
    query_b: str
    shared: int
    total_a: int
    total_b: int


def query_overlaps(index: Sequence[RecordRef]) -> List[QueryOverlap]:
    """
    Benign redundancy: how many of the same source records two queries both hit.

    Deliberately separate from collision reporting - shared records are expected
    when several queries search the same genome, and conflating the two would
    bury the cases that actually lose data.
    """
    hits: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
    for ref in index:
        hits[(ref.db, ref.query)].add(ref.source_id)

    dbs = sorted({db for db, _ in hits})
    overlaps: List[QueryOverlap] = []
    for db in dbs:
        queries = sorted(q for d, q in hits if d == db)
        for i, qa in enumerate(queries):
            for qb in queries[i + 1:]:
                a, b = hits[(db, qa)], hits[(db, qb)]
                shared = len(a & b)
                if shared:
                    overlaps.append(QueryOverlap(db=db, query_a=qa, query_b=qb,
                                                 shared=shared, total_a=len(a), total_b=len(b)))
    return overlaps


# -----------------------
# Console rendering
# -----------------------

def format_collision_lines(collisions: Sequence[Collision], limit: int = 10) -> List[str]:
    """
    Compact preview of what deduplication would do. Capped, because a gene-symbol
    header rule on a well-annotated genome routinely produces dozens of these and
    a wall of text is not a review.
    """
    lines: List[str] = []
    for col in collisions[:limit]:
        winner = col.winner
        lines.append(f"    {col.token}  <- keep {winner.source_id} ({winner.aa_len} aa)")
        for ref in col.losers:
            lines.append(f"        drop {ref.source_id} ({ref.aa_len} aa)")
    remaining = len(collisions) - limit
    if remaining > 0:
        lines.append(f"    ... and {remaining} more (full list in the de-duplication log)")
    return lines


def summarize(collisions: Sequence[Collision], res: Resolution) -> Dict[str, int]:
    """Counts for the end-of-run report."""
    counts = {ISOFORM: 0, CROSS_QUERY: 0, CROSS_DATABASE: 0,
              "declined": 0, "records_dropped": len(res.dropped)}
    for col in collisions:
        counts[col.kind] += 1
        if col.db is not None and (col.db, col.token) in res.declined_tokens:
            counts["declined"] += 1
    return counts


# -----------------------
# TSV log
# -----------------------

LOG_COLUMNS = [
    "stage", "query", "database", "identifier", "action",
    "source_id", "aa_len", "original_header", "reason",
]


def build_log_rows(index: Sequence[RecordRef], res: Resolution) -> List[Dict[str, object]]:
    """One row per affected record: what it was called, what happened, and why."""
    rows: List[Dict[str, object]] = []

    # A source record can be hit by several queries; name all of them rather
    # than whichever happened to be seen first.
    finders: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
    for ref in index:
        finders[(ref.db, ref.source_id)].add(ref.query)

    def row(stage, ref, action, identifier, reason, query=None):
        return {
            "stage": stage,
            "query": query or ",".join(sorted(finders[(ref.db, ref.source_id)])),
            "database": ref.db,
            "identifier": identifier,
            "action": action,
            "source_id": ref.source_id,
            "aa_len": ref.aa_len,
            "original_header": ref.description,
            "reason": reason,
        }

    for col in res.collisions:
        stage = STAGE_BY_KIND[col.kind]

        if col.kind == CROSS_DATABASE:
            for ref in col.members:
                # Members dropped by isoform resolution are logged under that
                # collision instead; only survivors get tagged here.
                final = res.final_ids.get((ref.db, ref.source_id))
                if final is None:
                    continue
                others = ", ".join(sorted({m.db for m in col.members if m.db != ref.db}))
                rows.append(row(stage, ref, "renamed", final,
                                f"identifier '{col.token}' also used by {others}; "
                                f"genome tag appended so no source is dropped"))
            continue

        declined = (col.db, col.token) in res.declined_tokens
        if declined:
            for ref in col.members:
                final = res.final_ids.get((ref.db, ref.source_id), col.token)
                rows.append(row(stage, ref, "renamed", final,
                                f"user declined deduplication of '{col.token}'; "
                                f"all {len(col.members)} records retained with rank suffixes"))
            continue

        winner = col.winner
        # Both rows carry the identifier the locus ends up with, so grepping the
        # final label finds the record that was dropped to make room for it.
        final = res.final_ids.get((winner.db, winner.source_id), col.token)
        rows.append(row(stage, winner, "kept", final,
                        f"longest amino-acid sequence ({winner.aa_len} aa) of "
                        f"{len(col.members)} records parsing to '{col.token}'"))
        for ref in col.losers:
            rows.append(row(stage, ref, "dropped", final,
                            f"{ref.aa_len} aa; superseded by {winner.source_id} "
                            f"({winner.aa_len} aa) for identifier '{col.token}'"))

    # Benign cross-query overlap: the identical source record hit by several
    # queries. Collapsed at merge time, recorded here so the ledger balances.
    seen: Dict[Tuple[str, str], str] = {}
    for ref in sorted(index, key=lambda r: (r.db, r.source_id, r.query)):
        key = (ref.db, ref.source_id)
        if key in res.dropped:
            continue
        if key in seen:
            rows.append(row("merge_overlap", ref, "dropped",
                            res.final_ids.get(key, ref.token),
                            f"identical record already contributed by query "
                            f"'{seen[key]}' (expected overlap, no data lost)",
                            query=ref.query))
        else:
            seen[key] = ref.query

    return rows


def write_log(
    path: Path,
    rows: Sequence[Dict[str, object]],
    *,
    entry: str,
    blast_type: str,
    hdr_rules_by_db: Dict[str, Tuple[str, str]],
    mode: str,
) -> Path:
    """Write the de-duplication ledger as TSV with a self-describing header."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        fh.write("# blast-align-tree de-duplication log\n")
        fh.write(f"# entry: {entry}\n")
        fh.write(f"# generated: {datetime.now().isoformat(timespec='seconds')}\n")
        fh.write(f"# blast_type: {blast_type}\n")
        fh.write(f"# duplicates mode: {mode}\n")
        fh.write(f"# policy: isoform/duplicated locus -> longest amino-acid sequence retained "
                 f"(ties: lexicographically smallest source ID); "
                 f"identifier shared across databases -> genome tag appended, nothing dropped\n")
        for db, (hdr, sfx) in sorted(hdr_rules_by_db.items()):
            rule = f"-hdr {hdr!r}" + (f" -hdr_sfx {sfx!r}" if sfx else "")
            fh.write(f"# header rule: {db}\t{rule}\n")
        writer = csv.DictWriter(fh, fieldnames=LOG_COLUMNS, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    return path

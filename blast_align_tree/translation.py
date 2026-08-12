#!/usr/bin/env python3
"""
Explicit, documented handling of stop codons when nucleotide hits are translated.

A tblastn run retrieves whole nucleotide records from the BLAST database and
translates them. Those records may contain in-frame stop codons - in a
pseudogene, in a frameshifted locus, or simply as the terminal stop of a normal
CDS. What the pipeline does with such a stop is a scientific decision, not an
implementation detail, so it is stated here and recorded per sequence.

Three policies are available:

  * **truncate** (default) - translation ends at the first in-frame stop. The
    reported protein is the truncated product the locus actually encodes, and
    every retained residue stays in exact register with the nucleotide record
    (residue *p* spans nucleotides ``3p-2..3p``).
  * **readthrough** - every codon is retained in register and each in-frame stop
    is written as ``X`` (unknown), so no residue is invented and downstream
    domain content is preserved for domain-architecture work.
  * **excise** - the v1.0 behaviour, retained only so published runs can be
    reproduced. Each stop codon is deleted and the downstream sequence joined to
    the upstream sequence. This yields a protein that the input genome does not
    encode, and it shifts every downstream residue by one position per excised
    stop, so feature coordinates drift. Not recommended for new work.

**Reading frame and strand are never inferred.** Translation always uses forward
frame +1 of the record exactly as stored in the database, which is correct for a
CDS database and wrong for a transcript database whose entries carry a 5' UTR.
Rather than silently re-framing, this module translates all six frames purely for
diagnosis and flags any record that reads further in another frame, so a
frame problem surfaces as a warning instead of as a quietly wrong protein.

Everything here is pure functions over strings and dataclasses, so the policy can
be tested without running BLAST.
"""

from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

try:
    from Bio.Seq import Seq
except ImportError:  # pragma: no cover - cli.py reports the missing dependency
    Seq = None


# -----------------------
# Policies
# -----------------------

TRUNCATE = "truncate"
READTHROUGH = "readthrough"
EXCISE = "excise"
POLICIES = (TRUNCATE, READTHROUGH, EXCISE)
DEFAULT_POLICY = TRUNCATE

POLICY_DESCRIPTION = {
    TRUNCATE: "translation ends at the first in-frame stop; the truncated protein is reported",
    READTHROUGH: "all codons retained in register; each in-frame stop written as 'X'",
    EXCISE: "legacy v1.0 behaviour: in-frame stop codons deleted and the flanking sequence joined",
}

# The frame that is actually translated. Never changes, never inferred.
USED_FRAME = "+1"
STRAND_NOTE = "+ (record as stored in the database; not inferred)"

# All six frames are translated for diagnosis only.
FRAMES = ("+1", "+2", "+3", "-1", "-2", "-3")

# Evidence threshold for saying "another frame reads this record better".
#
# In a wrong reading frame, stop codons occur at roughly 3/64 per codon, so
# stop-free runs average about 21 codons and long runs are rare. An alternative
# frame therefore only counts as evidence when it reads well past that length
# *and* clearly further than frame +1 does. Comparing longest stop-free runs -
# rather than requiring the alternative frame to be entirely stop-free - is what
# makes this work on real data: a genuine 5' UTR carries its own stop codons in
# the coding frame, so "no stops at all" would never fire on the very case this
# check exists to catch.
MIN_ALT_ORF_CODONS = 50
ALT_ORF_RATIO = 2.0


# -----------------------
# Flags
# -----------------------

INTERNAL_STOP = "internal_stop"
NO_START_CODON = "no_start_codon"
NO_TERMINAL_STOP = "no_terminal_stop"
LEN_NOT_MULTIPLE_OF_3 = "len_not_multiple_of_3"
FRAME_MISMATCH = "frame_mismatch"          # emitted as "frame_mismatch:+2"
POSSIBLE_UTR = "possible_utr"
POSSIBLE_WRONG_STRAND = "possible_wrong_strand"
EMPTY_TRANSLATION = "empty_translation"
ORF_MOSTLY_DOWNSTREAM = "orf_mostly_downstream"

FLAG_DESCRIPTION = {
    INTERNAL_STOP: "at least one in-frame stop codon before the end of the record",
    NO_START_CODON: "record does not begin with ATG",
    NO_TERMINAL_STOP: "record does not end with an in-frame stop codon",
    LEN_NOT_MULTIPLE_OF_3: "record length is not a multiple of 3; trailing bases were not translated",
    FRAME_MISMATCH: "another reading frame reads further without an internal stop",
    POSSIBLE_UTR: "a different forward frame reads further: the record may carry a 5' UTR",
    POSSIBLE_WRONG_STRAND: "a reverse frame reads further: the record may be stored on the opposite strand",
    EMPTY_TRANSLATION: "no complete codon could be translated",
    ORF_MOSTLY_DOWNSTREAM: ("more coding potential lies after the first stop than before it; "
                            "the truncated protein is unlikely to be the real product"),
}

# Flags that mean "the reported protein may not be the real reading of this locus".
FRAME_FLAGS = (POSSIBLE_UTR, POSSIBLE_WRONG_STRAND)


# -----------------------
# Sequence helpers
# -----------------------

def _frame_nt(seq: str, frame: str) -> str:
    """Nucleotides of one reading frame, trimmed to whole codons."""
    s = seq
    if frame.startswith("-"):
        s = str(Seq(s).reverse_complement())
    s = s[int(frame[1:]) - 1:]
    return s[: len(s) - (len(s) % 3)]


def _translate(nt: str) -> str:
    """Translate whole codons, keeping stops as '*'."""
    if not nt:
        return ""
    return str(Seq(nt).translate())


def _stop_positions(aa: str) -> List[int]:
    return [i for i, ch in enumerate(aa) if ch == "*"]


def _internal_stop_positions(aa: str) -> List[int]:
    """Stop positions excluding a stop in the final codon (a normal terminal stop)."""
    stops = _stop_positions(aa)
    if stops and stops[-1] == len(aa) - 1:
        stops = stops[:-1]
    return stops


def _longest_orf(aa: str) -> int:
    """Longest run of residues containing no stop."""
    return max((len(run) for run in aa.split("*")), default=0)


def excise_stop_codons(seq: str) -> str:
    """
    Legacy v1.0 excision, preserved verbatim as the single definition of that
    behaviour: delete every in-frame stop codon and join the flanking sequence.
    """
    stops = {"TAG", "TGA", "TAA"}
    bases = list(seq)
    i = 0
    while i + 2 < len(bases):
        if "".join(bases[i:i + 3]).upper() in stops:
            del bases[i:i + 3]
            # do not advance i; the next codon now sits at the same index
        else:
            i += 3
    return "".join(bases)


# -----------------------
# Data model
# -----------------------

@dataclass(frozen=True)
class TranslationStat:
    """What happened when one nucleotide record was translated."""
    source_id: str
    description: str
    nt_len: int
    start_codon: bool
    terminal_stop: bool
    n_internal_stops: int
    first_stop_aa_pos: Optional[int]        # 1-based, of the first *internal* stop
    n_codons: int
    aa_len_reported: int                    # length written under the policy in force
    aa_len_ranking: int                     # residues excluding stops; policy-independent
    aa_after_first_stop: int                # codons remaining after the first internal stop
    longest_downstream_orf_aa: int          # longest stop-free run after that stop
    best_frame: str
    flags: Tuple[str, ...]
    policy: str

    @property
    def has_internal_stop(self) -> bool:
        return self.n_internal_stops > 0

    @property
    def frame_suspect(self) -> bool:
        return any(f in FRAME_FLAGS for f in self.flags)

    def as_row(self) -> Dict[str, object]:
        return {
            "source_id": self.source_id,
            "original_header": self.description,
            "nt_len": self.nt_len,
            "frame": USED_FRAME,
            "strand": STRAND_NOTE,
            "start_codon": _yn(self.start_codon),
            "terminal_stop": _yn(self.terminal_stop),
            "n_internal_stops": self.n_internal_stops,
            "first_stop_aa_pos": "" if self.first_stop_aa_pos is None else self.first_stop_aa_pos,
            "n_codons": self.n_codons,
            "aa_len_reported": self.aa_len_reported,
            "aa_len_ranking": self.aa_len_ranking,
            "aa_after_first_stop": self.aa_after_first_stop,
            "longest_downstream_orf_aa": self.longest_downstream_orf_aa,
            "best_frame": self.best_frame,
            "flags": ";".join(self.flags),
            "policy": self.policy,
        }


def _yn(value: bool) -> str:
    return "yes" if value else "no"


# -----------------------
# Translation
# -----------------------

def translate_sequence(seq: str, policy: str = DEFAULT_POLICY) -> str:
    """Amino-acid sequence to report for one nucleotide record under `policy`."""
    if policy not in POLICIES:
        raise ValueError(f"unknown internal-stop policy: {policy!r}")

    seq = str(seq).upper()

    if policy == EXCISE:
        return _translate(_frame_nt(excise_stop_codons(seq), USED_FRAME))

    aa = _translate(_frame_nt(seq, USED_FRAME))

    if policy == TRUNCATE:
        # Ends at the first stop of any kind. For a well-formed CDS that is the
        # terminal stop, so a clean record is unaffected.
        return aa.split("*", 1)[0]

    # READTHROUGH: drop a terminal stop, mark internal ones as unknown residues.
    if aa.endswith("*"):
        aa = aa[:-1]
    return aa.replace("*", "X")


def analyze(seq: str, *, source_id: str = "", description: str = "",
            policy: str = DEFAULT_POLICY) -> TranslationStat:
    """Translate one record under `policy` and record what the stops implied."""
    seq = str(seq).upper()
    aa_by_frame = {f: _translate(_frame_nt(seq, f)) for f in FRAMES}
    aa = aa_by_frame[USED_FRAME]

    internal = _internal_stop_positions(aa)
    n_codons = len(aa)
    total_stops = len(_stop_positions(aa))

    if internal:
        first = internal[0]
        first_stop_aa_pos: Optional[int] = first + 1          # 1-based
        aa_after_first_stop = n_codons - (first + 1)
        longest_downstream = _longest_orf(aa[first + 1:])
    else:
        first_stop_aa_pos = None
        aa_after_first_stop = 0
        longest_downstream = 0

    best_frame = min(FRAMES, key=lambda f: (-_longest_orf(aa_by_frame[f]), FRAMES.index(f)))

    start_codon = seq[:3] == "ATG"
    terminal_stop = aa.endswith("*")

    # A record that starts with ATG, ends in a stop, and carries no internal stop
    # is a CDS in frame +1; there is nothing to second-guess. Short-circuiting
    # here is what keeps the frame flags informative - a well-formed CDS database
    # produces none of them at all.
    frame_1_is_clean_cds = start_codon and terminal_stop and not internal

    # Otherwise, look for a frame that reads substantially further than frame +1
    # (see MIN_ALT_ORF_CODONS / ALT_ORF_RATIO for why those bounds).
    if frame_1_is_clean_cds or not aa:
        alt_frame = None
    else:
        used_orf = _longest_orf(aa)
        cleaner = [
            f for f in FRAMES
            if f != USED_FRAME
            and _longest_orf(aa_by_frame[f]) >= MIN_ALT_ORF_CODONS
            and _longest_orf(aa_by_frame[f]) >= ALT_ORF_RATIO * used_orf
        ]
        alt_frame = (min(cleaner, key=lambda f: (-_longest_orf(aa_by_frame[f]), FRAMES.index(f)))
                     if cleaner else None)

    flags: List[str] = []
    if not aa:
        # Nothing was translated; the other observations would only be noise.
        flags.append(EMPTY_TRANSLATION)
        if len(seq) % 3:
            flags.append(LEN_NOT_MULTIPLE_OF_3)
    else:
        if internal:
            flags.append(INTERNAL_STOP)
        if not start_codon:
            flags.append(NO_START_CODON)
        if not terminal_stop:
            flags.append(NO_TERMINAL_STOP)
        if len(seq) % 3:
            flags.append(LEN_NOT_MULTIPLE_OF_3)
        if alt_frame is not None:
            flags.append(f"{FRAME_MISMATCH}:{alt_frame}")
            if alt_frame.startswith("-"):
                flags.append(POSSIBLE_WRONG_STRAND)
            else:
                flags.append(POSSIBLE_UTR)
        # Independent of frame: an in-frame 5' UTR, or any stop very early in the
        # record, leaves a short truncated protein with the bulk of the coding
        # potential still ahead of it. That reads as a normal frame +1 record to
        # the check above, so compare the two lengths directly.
        if internal and longest_downstream > first_stop_aa_pos - 1:
            flags.append(ORF_MOSTLY_DOWNSTREAM)

    return TranslationStat(
        source_id=source_id,
        description=description or source_id,
        nt_len=len(seq),
        start_codon=start_codon,
        terminal_stop=terminal_stop,
        n_internal_stops=len(internal),
        first_stop_aa_pos=first_stop_aa_pos,
        n_codons=n_codons,
        aa_len_reported=len(translate_sequence(seq, policy)),
        # Residues excluding every stop character. Identical under all three
        # policies, so isoform ranking cannot shift when the policy changes.
        aa_len_ranking=n_codons - total_stops,
        aa_after_first_stop=aa_after_first_stop,
        longest_downstream_orf_aa=longest_downstream,
        best_frame=best_frame,
        flags=tuple(flags),
        policy=policy,
    )


def translate_record(seq: str, *, source_id: str = "", description: str = "",
                     policy: str = DEFAULT_POLICY) -> Tuple[str, TranslationStat]:
    """Amino-acid sequence plus its stat, for the common case of needing both."""
    return translate_sequence(seq, policy), analyze(
        seq, source_id=source_id, description=description, policy=policy)


# -----------------------
# Per-(query, database) sidecars
# -----------------------

STAT_COLUMNS = [
    "source_id", "nt_len", "frame", "strand", "start_codon", "terminal_stop",
    "n_internal_stops", "first_stop_aa_pos", "n_codons", "aa_len_reported",
    "aa_len_ranking", "aa_after_first_stop", "longest_downstream_orf_aa",
    "best_frame", "flags", "policy", "original_header",
]
SIDECAR_COLUMNS = ["database", "query"] + STAT_COLUMNS
LOG_COLUMNS = ["identifier", "source_id", "database", "queries", "action"] + \
              [c for c in STAT_COLUMNS if c != "source_id"]

SIDECAR_SUFFIX = ".stops.tsv"


def write_sidecar(path: Path, stats: Sequence[TranslationStat], *,
                  database: str, query: str) -> Path:
    """
    One stat file per (query, database), written while that job's translation is
    still in hand. Sidecars rather than shared state because translation runs
    inside a thread pool; they are merged once every job has finished.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=SIDECAR_COLUMNS, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for stat in stats:
            row = stat.as_row()
            row.update({"database": database, "query": query})
            writer.writerow(row)
    return path


def read_sidecars(paths: Iterable[Path]) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for path in paths:
        if not Path(path).is_file():
            continue
        with Path(path).open(encoding="utf-8", newline="") as fh:
            rows.extend(csv.DictReader(fh, delimiter="\t"))
    return rows


def build_log_rows(
    sidecar_rows: Sequence[Dict[str, str]],
    final_ids: Dict[Tuple[str, str], str],
    dropped: Iterable[Tuple[str, str]],
) -> List[Dict[str, object]]:
    """
    Collapse per-query sidecar rows to one row per source record, labelled with
    the identifier that record ends up carrying in the tree.

    A record hit by several queries is translated identically every time, so the
    queries are joined rather than duplicated. Records that identifier
    reconciliation dropped keep a row, marked as such, so this report and the
    de-duplication log describe the same set of sequences.
    """
    dropped_set = set(dropped)
    merged: Dict[Tuple[str, str], Dict[str, object]] = {}
    queries: Dict[Tuple[str, str], set] = {}

    for row in sidecar_rows:
        key = (row["database"], row["source_id"])
        queries.setdefault(key, set()).add(row["query"])
        if key not in merged:
            merged[key] = {c: row.get(c, "") for c in STAT_COLUMNS}
            merged[key]["database"] = row["database"]

    out: List[Dict[str, object]] = []
    for key in sorted(merged):
        db, source_id = key
        row = dict(merged[key])
        row["queries"] = ",".join(sorted(queries[key]))
        row["identifier"] = final_ids.get(key, "")
        row["action"] = "dropped" if key in dropped_set else "kept"
        if not row["identifier"]:
            # Dropped by collision resolution: no tree label of its own.
            row["identifier"] = "-"
        out.append(row)
    return out


# -----------------------
# TSV report
# -----------------------

def write_log(path: Path, rows: Sequence[Dict[str, object]], *, entry: str,
              policy: str, blast_type: str) -> Path:
    """Write the per-sequence translation report with a self-describing header."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        fh.write("# blast-align-tree translation report\n")
        fh.write(f"# entry: {entry}\n")
        fh.write(f"# generated: {datetime.now().isoformat(timespec='seconds')}\n")
        fh.write(f"# blast_type: {blast_type}\n")
        fh.write(f"# internal stop policy: {policy} - {POLICY_DESCRIPTION[policy]}\n")
        fh.write(f"# reading frame: {USED_FRAME}, always; never inferred from the data\n")
        fh.write(f"# strand: {STRAND_NOTE}\n")
        fh.write("# all six frames are translated for diagnosis only and reported as best_frame\n")
        fh.write("# aa_len_ranking excludes stop characters and is identical under every "
                 "policy, so isoform ranking does not change with the policy\n")
        fh.write("# aa_after_first_stop and longest_downstream_orf_aa measure how much "
                 "reading frame remained after the first internal stop\n")
        for flag in (INTERNAL_STOP, NO_START_CODON, NO_TERMINAL_STOP,
                     LEN_NOT_MULTIPLE_OF_3, FRAME_MISMATCH, POSSIBLE_UTR,
                     POSSIBLE_WRONG_STRAND, ORF_MOSTLY_DOWNSTREAM, EMPTY_TRANSLATION):
            fh.write(f"# flag {flag}: {FLAG_DESCRIPTION[flag]}\n")
        writer = csv.DictWriter(fh, fieldnames=LOG_COLUMNS, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


# -----------------------
# Console summary
# -----------------------

def summarize(rows: Sequence[Dict[str, object]]) -> Dict[str, object]:
    """Counts for the end-of-run report."""
    flag_counts: Counter = Counter()
    by_db: Counter = Counter()
    flagged_by_db: Counter = Counter()

    for row in rows:
        db = str(row.get("database", ""))
        by_db[db] += 1
        flags = [f for f in str(row.get("flags", "")).split(";") if f]
        for flag in flags:
            flag_counts[flag.split(":", 1)[0]] += 1
        if any(f in FRAME_FLAGS for f in flags):
            flagged_by_db[db] += 1

    return {
        "total": len(rows),
        "flag_counts": flag_counts,
        "by_db": by_db,
        "frame_suspect_by_db": flagged_by_db,
    }


def print_report(rows: Sequence[Dict[str, object]], policy: str,
                 log_path: Optional[Path]) -> None:
    """Print what the stop policy did, loudly enough that a frame problem is seen."""
    stats = summarize(rows)
    total = stats["total"]
    counts: Counter = stats["flag_counts"]

    print(f"\n  [translation] {total} sequences translated in frame {USED_FRAME} "
          f"(strand: {STRAND_NOTE})")
    print(f"    Internal stop policy: {policy} - {POLICY_DESCRIPTION[policy]}")

    if not total:
        return

    with_stops = counts.get(INTERNAL_STOP, 0)
    if with_stops:
        verb = {TRUNCATE: "truncated at the first stop",
                READTHROUGH: "stops written as X",
                EXCISE: "stops excised"}[policy]
        print(f"    {with_stops} of {total} sequences contain an internal stop -> {verb}")
    else:
        print(f"    No internal stop codons found")

    if counts.get(ORF_MOSTLY_DOWNSTREAM):
        print(f"    {counts[ORF_MOSTLY_DOWNSTREAM]} flagged {ORF_MOSTLY_DOWNSTREAM}: more coding "
              f"potential after the first stop than before it")

    for flag in (LEN_NOT_MULTIPLE_OF_3, NO_TERMINAL_STOP, NO_START_CODON, EMPTY_TRANSLATION):
        if counts.get(flag):
            print(f"    {counts[flag]} flagged {flag}")

    suspect: Counter = stats["frame_suspect_by_db"]
    if suspect:
        print(f"\n    !! {sum(suspect.values())} sequences read further in another frame.")
        print(f"       Frame {USED_FRAME} is used regardless - these were NOT re-framed.")
        if counts.get(POSSIBLE_UTR):
            print(f"       {counts[POSSIBLE_UTR]} flagged {POSSIBLE_UTR}: the database may hold "
                  f"transcripts with 5' UTRs rather than CDS.")
        if counts.get(POSSIBLE_WRONG_STRAND):
            print(f"       {counts[POSSIBLE_WRONG_STRAND]} flagged {POSSIBLE_WRONG_STRAND}.")
        for db, n in sorted(suspect.items(), key=lambda kv: (-kv[1], kv[0])):
            print(f"       {n:>5}  {db}")
        print(f"       Check these databases before interpreting the affected proteins.")

    if log_path is not None:
        print(f"    Per-sequence report: {log_path}")

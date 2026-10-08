"""Tests for the homology report.

Instead of an e-value cutoff, each run reports the weakest hit every search
kept and whether the -n cap truncated that search. These cover choosing the
worst hit, the cap flag, and the written report, without a BLAST run.
"""

from blast_align_tree import cli
from blast_align_tree import genome_selector as gs
from blast_align_tree import homology


# sseqid pident ppos length qcovhsp evalue bitscore, as BLAST writes them.
HITS = (
    "AT5G45260.1\t100.000\t100.00\t1288\t100\t0.0\t2650\n"
    "AT4G12020.2\t35.200\t52.10\t610\t48\t3.21e-40\t160\n"
    "AT2G17060.1\t28.900\t46.30\t402\t31\t1.10e-08\t58.2\n"
)


def _hits(tmp_path, text=HITS):
    fp = tmp_path / ("q.db.seq.tblastn" + homology.HITS_SUFFIX)
    fp.write_text(text, encoding="utf-8")
    return homology.read_hits(fp)


def test_worst_hit_is_highest_evalue(tmp_path):
    row = homology.build_row("Q", "TAIR10cds.fa", "25", _hits(tmp_path),
                             {("TAIR10cds.fa", "AT2G17060.1"): "AT2G17060"})
    assert row["worst_hit"] == "AT2G17060.1"
    assert row["worst_hit_identifier"] == "AT2G17060"
    assert row["worst_evalue"] == "1.10e-08"
    assert row["worst_pct_identity"] == "28.9"
    assert row["worst_pct_similarity"] == "46.3"
    assert row["worst_query_coverage"] == 31
    assert row["best_evalue"] == "0"


def test_evalue_ties_break_on_lowest_bitscore(tmp_path):
    hits = _hits(tmp_path, "A\t40\t60\t100\t50\t1e-5\t45\nB\t40\t60\t100\t50\t1e-5\t41\n")
    assert homology.build_row("Q", "db", "25", hits, {})["worst_hit"] == "B"


def test_limited_by_n_only_when_cap_reached(tmp_path):
    hits = _hits(tmp_path)
    assert homology.build_row("Q", "db", "3", hits, {})["limited_by_n"] == "yes"
    assert homology.build_row("Q", "db", "25", hits, {})["limited_by_n"] == "no"


def test_dropped_worst_hit_has_no_identifier(tmp_path):
    row = homology.build_row("Q", "db", "25", _hits(tmp_path), {})
    assert row["worst_hit_identifier"] == "-"


def test_search_without_hits(tmp_path):
    row = homology.build_row("Q", "db", "25", [], {})
    assert row["hits_returned"] == 0
    assert row["limited_by_n"] == "no"
    assert row["worst_evalue"] == "-"


def test_report_round_trips_to_selector_summary(tmp_path):
    rows = [homology.build_row("Q", "TAIR10cds.fa", "3", _hits(tmp_path), {}),
            homology.build_row("Q", "Vung469cds.fa", "25", [], {})]
    homology.write_log(tmp_path / cli.HOMOLOGY_REPORT_NAME, rows,
                       entry="Q", blast_type="tblastn")
    assert gs.HOMOLOGY_REPORT_NAME == cli.HOMOLOGY_REPORT_NAME
    assert gs.log_stats(tmp_path) == [
        "Homology:  2 searches — weakest hit e-value 1.10e-08, 28.9% identity; "
        "1 reached the -n cap"
    ]

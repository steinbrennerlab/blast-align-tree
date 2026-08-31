"""Tests for the run provenance the Recent Runs tab reads.

The pipeline records the command that produced each run; the selector reads
that file back, summarizes the run's parameters and outputs, and offers both a
re-run and a re-draw command. These cover the round trip and the summaries,
without a BLAST run or a display.
"""

import sys

import pytest

from blast_align_tree import cli
from blast_align_tree import genome_selector as gs

COMMAND = [
    "blast-align-tree",
    "-q", "AT1G73080.1", "AT4G08850.1",
    "-qdbs", "TAIR10cds.fa", "TAIR10cds.fa",
    "-dbs", "TAIR10cds.fa", "Vung469cds.fa",
    "-hdr", "gene:", "locus=",
    "-n", "15", "20",
    "--tree_builder", "RAxML",
    "--datasets", "datasets/Vu_L2counts.txt",
]


@pytest.fixture
def run_dir(tmp_path, monkeypatch):
    """A run directory holding the files the pipeline archives."""
    monkeypatch.setattr(sys, "argv", COMMAND)
    cli.write_run_command(tmp_path / cli.RUN_COMMAND_NAME, "AT1G73080.1", "tblastn")

    assets = tmp_path / gs.RUN_ASSETS_DIRNAME
    assets.mkdir()
    # Support values sit on internal nodes; only tips should be counted.
    (assets / gs.COMBINED_TREE_NAME).write_text(
        "((A:0.1,B:0.2)0.9:0.3,(C:0.1,D:0.2)0.8:0.1);", encoding="utf-8")
    (assets / gs.MAPPING_NAME).write_text(
        "taxa\tgenome\n"
        "A\tTAIR10cds.fa\nB\tTAIR10cds.fa\nC\tVung469cds.fa\nD\tVung469cds.fa\n",
        encoding="utf-8")

    (tmp_path / gs.DEDUP_LOG_NAME).write_text(
        "# blast-align-tree de-duplication log\n"
        "stage\tquery\tdatabase\tidentifier\taction\n"
        "within_query\tQ\tTAIR10cds.fa\tX1\tkept\n"
        "within_query\tQ\tTAIR10cds.fa\tX2\tdropped\n"
        "within_query\tQ\tTAIR10cds.fa\tX3\tdropped\n", encoding="utf-8")
    (tmp_path / gs.TRANSLATION_REPORT_NAME).write_text(
        "# blast-align-tree translation report\n"
        "identifier\tsource_id\tn_internal_stops\n"
        "A\tA\t0\nB\tB\t2\n", encoding="utf-8")

    for name in ("AT1G73080.1.pdf", "AT1G73080.1.pdf.MSA.pdf",
                 "AT1G73080.1_redraw.pdf", "AT1G73080.1_redraw_XII.pdf"):
        (tmp_path / name).write_bytes(b"%PDF-1.4\n")
    return tmp_path


def test_recorded_command_round_trips(run_dir):
    """What the pipeline wrote is what the selector hands back to be re-run."""
    assert gs.read_run_command(run_dir) == " ".join(COMMAND)


def test_command_file_keeps_its_provenance_header(run_dir):
    header = (run_dir / cli.RUN_COMMAND_NAME).read_text(encoding="utf-8").splitlines()
    assert header[0] == "# blast-align-tree run command"
    assert "# entry: AT1G73080.1" in header
    assert "# blast_type: tblastn" in header


def test_missing_command_file_reads_as_none(tmp_path):
    """Runs archived before command logging have nothing to reconstruct from."""
    assert gs.read_run_command(tmp_path) is None


def test_describe_command_pairs_databases_with_their_n(run_dir):
    lines = gs.describe_command(gs.read_run_command(run_dir))
    assert "Queries:   AT1G73080.1, AT4G08850.1" in lines
    assert "Databases: TAIR10cds.fa (-n 15), Vung469cds.fa (-n 20)" in lines


def test_describe_command_marks_untouched_options_as_defaults(run_dir):
    lines = gs.describe_command(gs.read_run_command(run_dir))
    options = next(ln for ln in lines if ln.startswith("Options:"))
    assert "RAxML" in options
    assert "clustalo (default)" in options
    assert "tblastn (default)" in options


def test_command_flags_keeps_quoted_values_intact():
    flags = gs.command_flags("blast-align-tree -dbs TAIR10cds.fa -hdr '>' 'gene:'")
    assert flags["-dbs"] == ["TAIR10cds.fa"]
    assert flags["-hdr"] == [">", "gene:"]


def test_tree_stats_counts_tips_not_support_values(run_dir):
    assert gs.tree_stats(run_dir) == "4 tips — TAIR10cds.fa 2, Vung469cds.fa 2"


def test_log_stats_summarizes_both_reports(run_dir):
    lines = gs.log_stats(run_dir)
    assert lines[0] == "De-dup:    3 entries — 2 dropped, 1 kept"
    assert lines[1] == "Translated: 2 sequences — 1 with internal stops"


def test_log_stats_is_quiet_when_a_run_has_no_logs(tmp_path):
    assert gs.log_stats(tmp_path) == []


def test_redraw_names_skips_the_main_render_and_its_companions(run_dir):
    assert gs.redraw_names(run_dir, "AT1G73080.1") == [
        "AT1G73080.1_redraw", "AT1G73080.1_redraw_XII"]


def test_redraw_command_matches_the_pipeline_hint():
    cmd = gs.redraw_command("AT1G73080.1", "20260830_0930")
    assert '-e AT1G73080.1 -b AT1G73080.1_redraw' in cmd
    assert '--subdir "runs/20260830_0930" -n <NODE>' in cmd
    assert "visualize_tree.r" in cmd


def test_redraw_command_carries_datasets_from_the_original_run(run_dir):
    datasets = gs.command_flags(gs.read_run_command(run_dir))["--datasets"][0]
    cmd = gs.redraw_command("AT1G73080.1", "20260830_0930", datasets)
    assert cmd.endswith('--datasets "datasets/Vu_L2counts.txt"')


def test_format_timestamp_passes_through_unrecognized_names():
    assert gs.format_timestamp("20260131_0814") == "2026-01-31 08:14"
    assert gs.format_timestamp("scratch") == "scratch"

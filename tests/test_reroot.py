"""Tests for -a/--reroot, which roots the tree drawn at the end of a run.

visualize_tree.r matches the tip label exactly, so a mistyped ID would kill a
whole pipeline at its last step. These cover the tip lookup that prevents that,
and the plumbing from the flag to the Rscript command.
"""

import pytest

from blast_align_tree import cli
from blast_align_tree import genome_selector as gs

NEWICK = "((AT2G38240:0.1,Vigun08g171800:0.2)0.9:0.3,(Phvul.008G214200.1:0.1,AT1G05010:0.2)0.8:0.1);"


@pytest.fixture
def entry_dir(tmp_path):
    """An entry directory holding the tree as it exists when the R script runs."""
    (tmp_path / "combinedtree.nwk").write_text(NEWICK, encoding="utf-8")
    return tmp_path


def test_read_tip_labels_skips_support_values(entry_dir):
    assert cli.read_tip_labels(entry_dir / "combinedtree.nwk") == [
        "AT2G38240", "Vigun08g171800", "Phvul.008G214200.1", "AT1G05010"]


def test_read_tip_labels_is_quiet_when_there_is_no_tree(tmp_path):
    assert cli.read_tip_labels(tmp_path / "combinedtree.nwk") == []


def test_resolve_reroot_accepts_a_tip_of_the_tree(entry_dir):
    assert cli.resolve_reroot("AT2G38240", entry_dir) == "AT2G38240"


def test_resolve_reroot_drops_an_id_that_is_not_a_tip(entry_dir, capsys):
    assert cli.resolve_reroot("AT5G45250", entry_dir) is None
    assert "not a tip of the tree" in capsys.readouterr().out


def test_resolve_reroot_names_the_isoform_that_hdr_parsing_left_behind(entry_dir, capsys):
    # -hdr gene: strips '.1' from Arabidopsis tips, so the ID the user knows
    # the gene by is one suffix longer than the tip label.
    assert cli.resolve_reroot("AT2G38240.1", entry_dir) is None
    assert "Did you mean: AT2G38240?" in capsys.readouterr().out


def test_resolve_reroot_defers_to_the_r_script_when_no_tree_is_readable(tmp_path):
    assert cli.resolve_reroot("AT2G38240", tmp_path) == "AT2G38240"


def _rscript_command(monkeypatch, workdir, reroot):
    calls = []
    monkeypatch.setattr(cli, "run", lambda cmd, cwd=None, capture=False: calls.append(cmd))
    cli.visualize_tree("ENTRY", ["ENTRY"], workdir, None, reroot)
    return calls[0]


def test_visualize_tree_passes_a_valid_reroot_to_the_r_script(tmp_path, monkeypatch):
    (tmp_path / "ENTRY").mkdir()
    (tmp_path / "ENTRY" / "combinedtree.nwk").write_text(NEWICK, encoding="utf-8")
    cmd = _rscript_command(monkeypatch, tmp_path, "AT2G38240")
    assert cmd[-2:] == ["--reroot", "AT2G38240"]


def test_visualize_tree_still_draws_the_tree_when_the_reroot_is_unknown(tmp_path, monkeypatch):
    (tmp_path / "ENTRY").mkdir()
    (tmp_path / "ENTRY" / "combinedtree.nwk").write_text(NEWICK, encoding="utf-8")
    cmd = _rscript_command(monkeypatch, tmp_path, "AT5G45250")
    assert "--reroot" not in cmd
    assert cmd[-2:] == ["--write", "ENTRY"]


def test_visualize_tree_omits_reroot_when_the_flag_is_absent(tmp_path, monkeypatch):
    (tmp_path / "ENTRY").mkdir()
    cmd = _rscript_command(monkeypatch, tmp_path, None)
    assert "--reroot" not in cmd


def test_describe_command_surfaces_the_reroot_of_a_past_run():
    command = ("blast-align-tree -q AT2G19590.1 -qdbs TAIR10cds.fa "
               "-dbs TAIR10cds.fa -n 15 -add AT2G38240 --reroot AT2G38240")
    extras = next(ln for ln in gs.describe_command(command) if ln.startswith("Extras:"))
    assert "added seqs: AT2G38240" in extras
    assert "reroot: AT2G38240" in extras

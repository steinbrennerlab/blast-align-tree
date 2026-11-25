"""Launcher command-line options win over the saved session and set Reset All."""

import argparse
import json
import tkinter as tk

import pytest

from blast_align_tree import genome_selector as gs


def test_header_option_keeps_equals_signs_in_the_token():
    args = gs.parse_args([
        "--header", "Pvul218cds.fa=locus=", "--header", r"legumes\Gmax.fa=gene:",
        "--blast-type", "blastp", "--clear-queries", "--default-n", "60"])
    assert args.header == [("Pvul218cds.fa", "locus="), (r"legumes\Gmax.fa", "gene:")]
    assert args.blast_type == "blastp"
    assert args.clear_queries
    assert args.default_n == 60


@pytest.mark.parametrize("value", ["Pvul218cds.fa", "=locus=", "Pvul218cds.fa="])
def test_malformed_header_option_is_rejected(value):
    with pytest.raises(argparse.ArgumentTypeError):
        gs.genome_header(value)


def test_defaults_leave_the_saved_session_alone():
    args = gs.parse_args([])
    assert args.header == [] and args.blast_type is None and not args.clear_queries


@pytest.fixture(scope="module")
def tk_root():
    try:
        root = tk.Tk()
    except tk.TclError as exc:
        pytest.skip(f"Tk display unavailable: {exc}")
    root.withdraw()
    yield root
    root.destroy()


@pytest.fixture
def project(tmp_path, monkeypatch):
    genomes = tmp_path / "genomes"
    (genomes / "legumes").mkdir(parents=True)
    (genomes / "TAIR10cds.fa").write_text(
        ">AT1G01010.1 gene:AT1G01010 symbol:NAC001\nATG\n", encoding="utf-8")
    (genomes / "Pvul218cds.fa").write_text(
        ">Phvul.001G000100.1 locus=Phvul.001G000100 ID=x\nATG\n", encoding="utf-8")
    (genomes / "legumes" / "Gmax.fa").write_text(
        ">Glyma.01G000100.1 locus=Glyma.01G000100 ID=y\nATG\n", encoding="utf-8")
    (genomes / "ITAG4.1_CDS.fa").write_text(">Solyc01g005000.3.1\nATG\n",
                                             encoding="utf-8")
    monkeypatch.setattr(gs, "PROJ_DIR", tmp_path)
    monkeypatch.setattr(gs, "GENOMES_DIR", genomes)
    monkeypatch.setattr(gs, "SETTINGS_FILE", tmp_path / "settings.json")
    gs._id_cache.clear()
    return tmp_path


def make_app(tk_root, argv):
    args = gs.parse_args(argv)
    window = tk.Toplevel(tk_root)
    window.withdraw()
    app = gs.GenomeSelectorApp(
        window, default_n=args.default_n, selected_genomes=args.select_genome,
        restore=not args.no_restore, headers=args.header,
        blast_type=args.blast_type, clear_queries=args.clear_queries)
    return app, window


LAUNCHER = [
    "--default-n", "60", "--blast-type", "tblastn", "--clear-queries",
    "--select-genome", "TAIR10cds.fa", "--select-genome", "Pvul218cds.fa",
    "--select-genome", "legumes/Gmax.fa",
    "--header", "TAIR10cds.fa=gene:", "--header", "Pvul218cds.fa=locus=",
    "--header", "legumes/Gmax.fa=locus=",
]
EXPECTED_ROWS = {
    "TAIR10cds.fa": {"selected": True, "hdr": "gene:", "sfx": "none", "n": "60", "q": ""},
    "Pvul218cds.fa": {"selected": True, "hdr": "locus=", "sfx": "none", "n": "60", "q": ""},
    "legumes/Gmax.fa": {"selected": True, "hdr": "locus=", "sfx": "none", "n": "60", "q": ""},
    "ITAG4.1_CDS.fa": {"selected": False, "hdr": "id", "sfx": "none", "n": "60", "q": ""},
}


def test_launcher_options_override_the_saved_session(project, tk_root):
    saved = {
        "version": gs.SETTINGS_VERSION,
        "globals": {"default_n": "40", "blast_type": "blastp", "aligner": "mafft"},
        "rows": {
            "ITAG4.1_CDS.fa": {"selected": True, "hdr": "id", "sfx": "none",
                               "n": "40", "q": "Solyc01g005000.3.1"},
            "TAIR10cds.fa": {"selected": True, "hdr": "id", "sfx": "none",
                             "n": "40", "q": "AT1G01010.1"},
            "Pvul218cds.fa": {"selected": False, "hdr": "ID=", "sfx": "none",
                              "n": "40", "q": ""},
        },
    }
    gs.SETTINGS_FILE.write_text(json.dumps(saved), encoding="utf-8")
    app, window = make_app(tk_root, LAUNCHER)
    try:
        assert app._snapshot_rows() == EXPECTED_ROWS
        assert app.default_n_var.get() == "60"
        assert app.blast_type_var.get() == "tblastn"
        # Options the launcher does not name still come from the saved session.
        assert app.aligner_var.get() == "mafft"
        assert app.output.get("1.0", "end").strip() == ""

        # Reset All returns to the launcher's state, not the saved session's.
        app.default_n_var.set("5")
        app._set_all_n()
        app.rows[0][4].set("typed query")
        app._apply_settings(app._launch_defaults)
        assert app._snapshot_rows() == EXPECTED_ROWS
        assert app.blast_type_var.get() == "tblastn"
    finally:
        window.destroy()


def test_launcher_options_apply_without_a_saved_session(project, tk_root):
    app, window = make_app(tk_root, LAUNCHER + ["--no-restore"])
    try:
        assert app._snapshot_rows() == EXPECTED_ROWS
        assert app._launch_defaults["rows"] == EXPECTED_ROWS
    finally:
        window.destroy()


def test_unusable_launcher_options_are_reported(project, tk_root):
    app, window = make_app(tk_root, [
        "--no-restore", "--header", "Pvul218cds.fa=gene:",
        "--header", "missing.fa=id", "--select-genome", "Gone.fa"])
    try:
        assert app._snapshot_rows()["Pvul218cds.fa"]["hdr"] == "id"
        output = app.output.get("1.0", "end")
        assert "Pvul218cds.fa: header 'gene:' is not in this file" in output
        assert "missing.fa: no such genome" in output
        assert "gone.fa: no such genome" in output
    finally:
        window.destroy()

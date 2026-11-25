"""Loading a run restores editable settings without carrying over another run."""

from pathlib import Path
import tkinter as tk

import pytest

from blast_align_tree import genome_selector as gs


COMMAND = (
    'blast-align-tree -q Q1 Q2 Q3 -qdbs a.fa a.fa nested/b.fa '
    '-dbs nested/b.fa a.fa -hdr "gene:" id -hdr_sfx .1 none -n 20 5 '
    '--aligner mafft --mafft_mode linsi --tree_builder RAxML '
    '--blast_type blastp --threads=3 -add OUT -add_db out.fa -a OUT '
    '-aa 10 200 --motif "M=A[ST]" --motif_syntax regex --motif_overlap '
    '--hmm profiles/kinase.hmm'
)


def test_loads_database_pairs_queries_and_advanced_options():
    settings, warnings = gs.settings_from_command(
        COMMAND, ["a.fa", "nested/b.fa", "unused.fa"])
    assert warnings == []
    assert settings["rows"]["a.fa"] == {
        "selected": True, "hdr": "id", "sfx": "none", "n": "5", "q": "Q1 Q2"}
    assert settings["rows"]["nested/b.fa"] == {
        "selected": True, "hdr": "gene:", "sfx": ".1", "n": "20", "q": "Q3"}
    assert not settings["rows"]["unused.fa"]["selected"]
    assert settings["rows"]["unused.fa"]["q"] == ""
    assert settings["globals"] == {
        "default_n": "15", "prepend": True, "aligner": "mafft",
        "mafft_mode": "linsi", "tree_builder": "RAxML", "blast_type": "blastp",
        "threads": "3", "add_seqs": "OUT", "add_dbs": "out.fa", "reroot": "OUT",
        "aa_start": "10", "aa_end": "200", "motif": "M=A[ST]",
        "motif_syntax": "regex", "motif_overlap": True,
        "hmm": "profiles/kinase.hmm", "adv_visible": True,
    }


def test_long_aliases_and_query_only_database():
    settings, warnings = gs.settings_from_command(
        'blast-align-tree --queries Q --query_databases query.fa '
        '--database hits.fa --header custom: --header_suffix .x --seqs 7 '
        '--add_seqs OUT --add_dbs out.fa --reroot OUT --slice 2 80',
        ["hits.fa", "query.fa"])
    assert not warnings
    assert settings["rows"]["query.fa"]["q"] == "Q"
    assert not settings["rows"]["query.fa"]["selected"]
    assert settings["rows"]["hits.fa"]["sfx"] == ".x"
    assert settings["globals"]["add_seqs"] == "OUT"
    assert settings["globals"]["aa_end"] == "80"


def test_missing_genomes_and_unsupported_options_are_reported():
    settings, warnings = gs.settings_from_command(
        COMMAND + ' --datasets data.tsv --duplicates auto', ["a.fa", "b.fa"])
    assert "Genome not available in the panel: nested/b.fa" in warnings
    assert not settings["rows"]["b.fa"]["selected"]
    assert any("--datasets" in warning and "--duplicates" in warning for warning in warnings)


def test_relative_path_separators_are_normalized():
    settings, warnings = gs.settings_from_command(
        COMMAND.replace("nested/b.fa", r".\nested\b.fa"), ["a.fa", "nested/b.fa"])
    assert not warnings
    assert settings["rows"]["nested/b.fa"]["q"] == "Q3"


def test_per_query_slices_are_reported_without_loading_an_incorrect_range():
    settings, warnings = gs.settings_from_command(
        COMMAND.replace("-aa 10 200", "-aa 10:200 none 2:40"), ["a.fa", "nested/b.fa"])
    assert settings["globals"]["aa_start"] == settings["globals"]["aa_end"] == ""
    assert any("slices" in warning for warning in warnings)


@pytest.mark.parametrize("command", [
    "", "blast-align-tree -q Q", COMMAND.replace("-n 20 5", "-n 20"),
    COMMAND.replace("-qdbs a.fa a.fa nested/b.fa", "-qdbs a.fa"),
    COMMAND + ' --hmm "unterminated',
])
def test_malformed_commands_do_not_produce_partial_settings(command):
    with pytest.raises(ValueError):
        gs.settings_from_command(command, ["a.fa", "nested/b.fa"])


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
def app(tmp_path, monkeypatch, tk_root):
    root = tk_root
    monkeypatch.setattr(gs, "PROJ_DIR", tmp_path)
    monkeypatch.setattr(gs, "SETTINGS_FILE", tmp_path / "settings.json")
    app = gs.GenomeSelectorApp.__new__(gs.GenomeSelectorApp)
    app.root = root
    app.default_n = 99
    app.row_tokens = {"a.fa": ["id"], "nested/b.fa": ["id"], "unused.fa": ["id"]}
    app.rows = [(tk.BooleanVar(root, True), tk.StringVar(root, "stale:"),
                 tk.StringVar(root, ".stale"), tk.StringVar(root, "99"),
                 tk.StringVar(root, "STALE"), Path(name), name)
                for name in app.row_tokens]
    for key, attr in gs.GLOBAL_SETTING_VARS.items():
        var_type = tk.BooleanVar if key in ("motif_overlap", "prepend") else tk.StringVar
        setattr(app, attr, var_type(root))
    app._adv_visible = tk.BooleanVar(root, False)
    monkeypatch.setattr(app, "_toggle_advanced",
                        lambda: app._adv_visible.set(not app._adv_visible.get()))
    app.selected_label = tk.Label(root)
    app._build_output_panel()
    yield app
    app.notebook.destroy()
    app.selected_label.destroy()


def select_run(app, tmp_path, command):
    path = tmp_path / "Q1" / "runs" / "20260922_1200"
    path.mkdir(parents=True, exist_ok=True)
    (path / gs.RUN_COMMAND_NAME).write_text(command, encoding="utf-8")
    app._refresh_runs()
    return path


def test_load_button_replaces_stale_fields_and_round_trips(app, tmp_path, monkeypatch):
    path = select_run(app, tmp_path, COMMAND)
    app.load_run_options_btn.invoke()
    assert app._snapshot_rows()["nested/b.fa"]["hdr"] == "gene:"
    assert app._snapshot_rows()["unused.fa"]["q"] == ""
    assert not app._snapshot_rows()["unused.fa"]["selected"]
    assert app._adv_visible.get()
    assert "Loaded options from Q1" in app.output.get("1.0", "end")
    assert app.notebook.index(app.notebook.select()) == 0
    monkeypatch.setattr(gs, "load_fasta_ids", lambda path: {"Q1", "Q2", "Q3"})
    monkeypatch.setattr(gs, "dbtype_for", lambda path: "prot")
    rebuilt = "blast-align-tree " + " ".join(app._build_args())
    expected, _ = gs.settings_from_command(COMMAND, list(app.row_tokens))
    assert gs.settings_from_command(rebuilt, list(app.row_tokens))[0] == expected
    assert (path / gs.RUN_COMMAND_NAME).read_text(encoding="utf-8") == COMMAND

    # A second run with omitted flags must reset advanced and global options.
    select_run(app, tmp_path, "blast-align-tree -q Q1 -qdbs a.fa -dbs a.fa -hdr id -n 8")
    app.load_run_options_btn.invoke()
    assert app.aligner_var.get() == "clustalo"
    assert app.tree_var.get() == "FastTree"
    assert app.blast_type_var.get() == "tblastn"
    assert app.add_seqs_var.get() == app.aa_start_var.get() == app.hmm_var.get() == ""
    assert not app.motif_overlap_var.get()
    assert not app._adv_visible.get()
    assert not app._snapshot_rows()["nested/b.fa"]["selected"]


def test_missing_command_disables_load_and_invalid_command_preserves_form(app, tmp_path):
    assert app.load_run_options_btn.instate(["disabled"])
    path = select_run(app, tmp_path, "")
    assert app.load_run_options_btn.instate(["disabled"])
    select_run(app, tmp_path, "blast-align-tree -q Q")
    before = app._snapshot()
    app.load_run_options_btn.invoke()
    assert app._snapshot() == before
    assert "Could not load run options" in app.output.get("1.0", "end")
    (path / gs.RUN_COMMAND_NAME).unlink()
    app._refresh_runs()
    assert app.load_run_options_btn.instate(["disabled"])

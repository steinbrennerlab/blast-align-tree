"""Run labels persist independently of pipeline outputs and ordering."""

import tkinter as tk
from pathlib import Path
from tkinter import ttk

import pytest

from blast_align_tree import genome_selector as gs


@pytest.fixture
def runs(tmp_path, monkeypatch):
    monkeypatch.setattr(gs, "PROJ_DIR", tmp_path)
    result = []
    for query, timestamp in [("zQuery", "20260922_1200"),
                              ("aQuery", "20260921_1200"),
                              ("aQuery", "20260920_1200")]:
        path = tmp_path / query / "runs" / timestamp
        path.mkdir(parents=True)
        result.append((query, timestamp, path))
    return result


def test_metadata_round_trip_is_specific_to_one_run(runs):
    path = runs[1][2]
    command = path / gs.RUN_COMMAND_NAME
    command.write_text("blast-align-tree -q aQuery", encoding="utf-8")
    metadata = {"name": "Receptor α", "starred": True}
    gs.write_run_metadata(path, metadata)
    assert gs.read_run_metadata(path) == metadata
    assert gs.read_run_metadata(runs[2][2]) == {"name": "", "starred": False}
    assert command.read_text(encoding="utf-8") == "blast-align-tree -q aQuery"
    gs.write_run_metadata(path, {"name": "", "starred": False})
    assert gs.read_run_metadata(path) == {"name": "", "starred": False}


@pytest.mark.parametrize("contents", ["{broken", "[]", "null",
                                         '{"name": 4, "starred": "false"}'])
def test_invalid_metadata_is_safe(tmp_path, contents):
    (tmp_path / gs.RUN_METADATA_NAME).write_text(contents, encoding="utf-8")
    assert gs.read_run_metadata(tmp_path) == {"name": "", "starred": False}


def test_failed_save_preserves_previous_metadata(tmp_path, monkeypatch):
    original = {"name": "Keep me", "starred": True}
    gs.write_run_metadata(tmp_path, original)

    def fail_replace(*args):
        raise PermissionError("Read-only run")

    monkeypatch.setattr(Path, "replace", fail_replace)
    with pytest.raises(PermissionError):
        gs.write_run_metadata(tmp_path, {"name": "New name", "starred": False})
    assert gs.read_run_metadata(tmp_path) == original
    assert list(tmp_path.iterdir()) == [tmp_path / gs.RUN_METADATA_NAME]


@pytest.mark.parametrize("sort_by, indices", [
    ("Newest first", [0, 1, 2]),
    ("Oldest first", [2, 1, 0]),
    ("Starred first", [2, 0, 1]),
    ("Query A–Z", [1, 2, 0]),
    ("Query Z–A", [0, 1, 2]),
    ("Name A–Z", [2, 1, 0]),
])
def test_sort_modes_only_prioritize_stars_when_requested(runs, sort_by, indices):
    metadata = {runs[1][2]: {"name": "Beta"},
                runs[2][2]: {"name": "Alpha", "starred": True}}
    assert gs.sort_runs(list(reversed(runs)), sort_by, metadata) == [
        runs[index] for index in indices]


def test_star_sort_includes_runs_older_than_the_previous_limit(tmp_path, monkeypatch):
    monkeypatch.setattr(gs, "PROJ_DIR", tmp_path)
    for day in range(1, 26):
        (tmp_path / "Query" / "runs" / f"202608{day:02d}_1200").mkdir(parents=True)
    runs = gs.scan_recent_runs(limit=None)
    assert len(runs) == 25
    oldest = runs[-1]
    gs.write_run_metadata(oldest[2], {"name": "Old favorite", "starred": True})
    metadata = {run[2]: gs.read_run_metadata(run[2]) for run in runs}
    assert gs.sort_runs(runs, "Starred first", metadata)[0] == oldest
    assert gs.sort_runs(runs, "Newest first", metadata)[-1] == oldest


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
def app(runs, tk_root):
    app = gs.GenomeSelectorApp.__new__(gs.GenomeSelectorApp)
    app.root = tk_root
    app.notebook = ttk.Notebook(tk_root)
    app._build_recent_runs_tab()
    tk_root.update_idletasks()
    yield app
    app.notebook.destroy()


def displayed_queries(app):
    labels = [widget for widget in app.runs_inner.winfo_children()
              if widget.grid_info().get("column") == 2
              and widget.grid_info().get("row", 0) > 0]
    return [widget.cget("text") for widget in labels]


def test_run_controls_persist_and_keep_selection(app, runs, monkeypatch):
    app.runs_inner.grid_slaves(row=3, column=0)[0].invoke()
    assert displayed_queries(app) == ["zQuery", "aQuery", "aQuery"]
    assert app._selected_run == runs[2]
    monkeypatch.setattr(gs.simpledialog, "askstring", lambda *a, **kw: "  Favorite α  ")
    app.runs_inner.grid_slaves(row=3, column=1)[0].invoke()
    assert gs.read_run_metadata(runs[2][2]) == {"name": "Favorite α", "starred": True}
    assert "Favorite α" in app.runs_info.get("1.0", "end")
    app.runs_sort_var.set("Starred first")
    app._refresh_runs()
    assert displayed_queries(app) == ["aQuery", "zQuery", "aQuery"]
    assert app._selected_run == runs[2]
    app.runs_sort_var.set("Query A–Z")
    app._refresh_runs()
    assert displayed_queries(app) == ["aQuery", "aQuery", "zQuery"]
    assert app._selected_run == runs[2]
    app._build_recent_runs_tab()  # A fresh tab reads the saved labels.
    assert app._run_metadata[runs[2][2]] == {"name": "Favorite α", "starred": True}


def test_cancel_and_clear_name_and_unstar(app, runs, monkeypatch):
    app._update_run_metadata(runs[1], name="Keep", starred=True)
    monkeypatch.setattr(gs.simpledialog, "askstring", lambda *a, **kw: None)
    app._rename_run(runs[1])
    assert gs.read_run_metadata(runs[1][2])["name"] == "Keep"
    monkeypatch.setattr(gs.simpledialog, "askstring", lambda *a, **kw: "   ")
    app._rename_run(runs[1])
    app._toggle_run_star(runs[1])
    assert gs.read_run_metadata(runs[1][2]) == {"name": "", "starred": False}


def test_save_error_is_reported_without_changing_ui(app, runs, monkeypatch):
    errors = []
    monkeypatch.setattr(gs.messagebox, "showerror", lambda *a, **kw: errors.append(a))

    def fail_save(*args):
        raise PermissionError("Read-only run")

    monkeypatch.setattr(gs, "write_run_metadata", fail_save)
    app._toggle_run_star(runs[0])
    assert errors and "Read-only run" in errors[0][1]
    assert not app._run_metadata[runs[0][2]]["starred"]

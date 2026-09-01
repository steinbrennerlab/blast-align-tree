#!/usr/bin/env python3
"""
BAT Genome Selector — tkinter GUI for building blast-align-tree commands.

Scans genomes/ for *.fa and *.fna files, auto-detects header tokens from the
first FASTA record, and generates copy-paste-ready -dbs / -hdr / -n arguments.
"""

import argparse
import json
import os
import platform
import shlex
import shutil
import random
import re
import subprocess
import tkinter as tk
from tkinter import messagebox, ttk
from pathlib import Path
from importlib.resources import files as _pkg_files


PROJ_DIR = Path.cwd()
# The tree-drawing R script, resolved the same way cli.py resolves it, so the
# re-draw hint below matches the one the pipeline prints when it finishes.
VISUALIZE_TREE_R = Path(str(_pkg_files("blast_align_tree") / "data" / "visualize_tree.r"))
# Files the pipeline leaves in an archived run. Names mirror cli.py; they are
# repeated rather than imported so the GUI stays free of the pipeline's
# Biopython/R dependencies.
RUN_COMMAND_NAME = "run_command.txt"
DEDUP_LOG_NAME = "deduplication_log.tsv"
TRANSLATION_REPORT_NAME = "translation_report.tsv"
RUN_ASSETS_DIRNAME = "genes_alignments_trees"
MAPPING_NAME = "merged_genome_mapping.txt"
COMBINED_TREE_NAME = "combinedtree.nwk"
# Pipeline defaults, so run details can say which options were left alone.
CLI_DEFAULTS = {"--blast_type": "tblastn", "--aligner": "clustalo",
                "--tree_builder": "FastTree"}
# The notebook only gets the height the panels above it leave over, so the
# Recent Runs tab asks for a modest share and scrolls within it.
RUNS_LIST_HEIGHT = 120
DETAILS_MAX_LINES = 12
GENOMES_DIR = PROJ_DIR / "genomes"
SETTINGS_FILE = PROJ_DIR / ".bat_selector_settings.json"
SETTINGS_VERSION = 1

# Saved-settings key -> attribute holding the tk variable. Per-row state is
# handled separately, keyed by each genome's path relative to genomes/.
GLOBAL_SETTING_VARS = {
    "default_n": "default_n_var",
    "aligner": "aligner_var",
    "mafft_mode": "mafft_mode_var",
    "tree_builder": "tree_var",
    "blast_type": "blast_type_var",
    "threads": "threads_var",
    "add_seqs": "add_seqs_var",
    "add_dbs": "add_dbs_var",
    "reroot": "reroot_var",
    "aa_start": "aa_start_var",
    "aa_end": "aa_end_var",
    "motif": "motif_var",
    "motif_syntax": "motif_syntax_var",
    "motif_overlap": "motif_overlap_var",
    "hmm": "hmm_var",
    "prepend": "prepend_var",
}
FASTA_EXTENSIONS = {".fa", ".faa", ".fas", ".fasta", ".fna"}
# BLAST index extensions to exclude
BLAST_INDEX_EXTS = {".ndb", ".nhr", ".nin", ".njs", ".nog", ".nos", ".not",
                    ".nsq", ".ntf", ".nto", ".pdb", ".phr", ".pin", ".pjs",
                    ".pog", ".pos", ".pot", ".psq", ".ptf", ".pto"}
MAX_LABEL_LEN = 40
DEFAULT_N = 15
MIN_N = 1
MAX_N = 999
DIVIDER_COL = 3
DIVIDER_COL_WIDTH = 12
# Cache of FASTA record IDs per genome file {filepath: set_of_ids}
_id_cache: dict[Path, set[str]] = {}


def positive_int(value: str) -> int:
    """Parse a positive integer for CLI options."""
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a positive integer") from exc
    if parsed < MIN_N:
        raise argparse.ArgumentTypeError("must be a positive integer")
    if parsed > MAX_N:
        raise argparse.ArgumentTypeError(f"must be {MAX_N} or less")
    return parsed


def has_blastdb(filepath: Path) -> bool:
    """Check whether BLAST database index files exist for a FASTA file."""
    return Path(str(filepath) + ".nhr").exists() or Path(str(filepath) + ".phr").exists()


_NUCL_CHARS = set("ATGCNUatgcnu")


def dbtype_for(filepath: Path) -> str:
    """Detect whether a FASTA file contains protein or nucleotide sequences."""
    if filepath.suffix == ".faa":
        return "prot"
    seq_chars: list[str] = []
    try:
        with open(filepath, encoding="utf-8", errors="replace") as f:
            in_seq = False
            for line in f:
                if line.startswith(">"):
                    if in_seq:
                        break
                    in_seq = True
                    continue
                if in_seq:
                    seq_chars.extend(line.strip())
    except OSError:
        return "nucl"
    if not seq_chars:
        return "nucl"
    nucl_frac = sum(1 for c in seq_chars if c in _NUCL_CHARS) / len(seq_chars)
    return "nucl" if nucl_frac > 0.8 else "prot"


def get_example_id(filepath: Path, target: int = 500) -> str:
    """Return the identifier from the target-th FASTA record (1-indexed)."""
    first_id = None
    count = 0
    try:
        with open(filepath, encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith(">"):
                    count += 1
                    record_id = line[1:].split()[0]
                    if count == 1:
                        first_id = record_id
                    if count == target:
                        return record_id
    except OSError:
        pass
    return first_id or ""


def get_example_header(filepath: Path, target: int = 500) -> str:
    """Return the full header line (without '>') from the target-th FASTA record."""
    first_hdr = None
    count = 0
    try:
        with open(filepath, encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith(">"):
                    count += 1
                    hdr = line[1:].strip()
                    if count == 1:
                        first_hdr = hdr
                    if count == target:
                        return hdr
    except OSError:
        pass
    return first_hdr or ""


def parse_header_token(description: str, headerword: str, suffix: str = "") -> str:
    """Preview the parsed name for a header line given a -hdr token and optional suffix."""
    fallback_id = description.split()[0] if description else ""
    if headerword == "id":
        token = fallback_id
    elif headerword not in description:
        token = fallback_id
    else:
        try:
            part = description.split(headerword, 1)[1]
            token = part.split(" ", 1)[0]
        except Exception:
            token = fallback_id
    if suffix and suffix != "none" and token.endswith(suffix):
        token = token[:-len(suffix)]
    return token


def detect_header_tokens(filepath: Path) -> list[str]:
    """Read the first FASTA header line and extract candidate -hdr tokens."""
    try:
        with open(filepath, encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith(">"):
                    header = line[1:].strip()
                    break
            else:
                return ["id"]
    except OSError:
        return ["id"]

    tokens = []
    seen = set()

    for m in re.finditer(r'\[(\w+)=', header):
        tok = "[" + m.group(1) + "="
        if tok not in seen:
            seen.add(tok)
            tokens.append(tok)

    stripped = re.sub(r'\[[^\]]*\]', '', header)

    for m in re.finditer(r'\b([A-Za-z_]\w*[=:])', stripped):
        tok = m.group(1)
        if tok not in seen:
            seen.add(tok)
            tokens.append(tok)

    if "id" not in seen:
        tokens.append("id")
    return tokens


def scan_genomes() -> list[Path]:
    """Return sorted list of FASTA files in genomes/ (including subfolders)."""
    if not GENOMES_DIR.is_dir():
        return []
    files = [
        p for p in GENOMES_DIR.rglob("*")
        if p.is_file() and p.suffix in FASTA_EXTENSIONS and p.suffix not in BLAST_INDEX_EXTS
    ]
    files.sort(key=lambda p: p.relative_to(GENOMES_DIR).as_posix().lower())
    return files


def load_fasta_ids(filepath: Path) -> set[str]:
    """Return the set of record IDs in a FASTA file. Results are cached."""
    if filepath in _id_cache:
        return _id_cache[filepath]
    ids: set[str] = set()
    try:
        with open(filepath, encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith(">"):
                    ids.add(line[1:].split()[0])
    except OSError:
        pass
    _id_cache[filepath] = ids
    return ids


def truncate_name(name: str, maxlen: int = MAX_LABEL_LEN) -> str:
    if len(name) <= maxlen:
        return name
    return name[: maxlen - 3] + "..."


def scan_recent_runs(limit: int = 20) -> list[tuple[str, str, Path]]:
    """Scan for pipeline run directories (ENTRY/runs/TIMESTAMP/).

    Returns a list of (entry, timestamp, full_path) sorted newest-first,
    capped at *limit*.
    """
    results = []
    for entry_dir in PROJ_DIR.iterdir():
        if not entry_dir.is_dir():
            continue
        runs_dir = entry_dir / "runs"
        if not runs_dir.is_dir():
            continue
        for run_dir in runs_dir.iterdir():
            if run_dir.is_dir():
                results.append((entry_dir.name, run_dir.name, run_dir))
    # Sort by timestamp descending (newest first)
    results.sort(key=lambda r: r[1], reverse=True)
    return results[:limit]


def redraw_command(entry: str, timestamp: str, datasets: str | None = None) -> str:
    """Return the Rscript command that re-draws the trees of an archived run.

    Mirrors the "To re-draw trees" hint the pipeline prints when it finishes.
    The entry doubles as the -b basename because cli.py takes both from the
    first query; <NODE> is left as a placeholder for the user.
    """
    cmd = (f'Rscript "{VISUALIZE_TREE_R}" -e {entry} -b {entry}_redraw '
           f'--subdir "runs/{timestamp}" -n <NODE>')
    if datasets:
        cmd += f' --datasets "{datasets}"'
    return cmd


def format_timestamp(timestamp: str) -> str:
    """20260131_0814 -> 2026-01-31 08:14; unrecognized names pass through."""
    if len(timestamp) >= 13 and "_" in timestamp:
        return (f"{timestamp[:4]}-{timestamp[4:6]}-{timestamp[6:8]}"
                f" {timestamp[9:11]}:{timestamp[11:13]}")
    return timestamp


def read_run_command(run_dir: Path) -> str | None:
    """Return the invocation recorded for a run, or None if it has none.

    Runs archived before the pipeline started recording commands have no such
    file, and nothing on disk records -n, the aligner, or -add, so there is
    nothing to fall back on.
    """
    try:
        text = (run_dir / RUN_COMMAND_NAME).read_text(encoding="utf-8",
                                                      errors="replace")
    except OSError:
        return None
    for line in text.splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            return line
    return None


def _is_flag(token: str) -> bool:
    return len(token) > 1 and token[0] == "-" and (token[1].isalpha() or token[1] == "-")


def command_flags(command: str) -> dict[str, list[str]]:
    """Group a recorded command into {flag: [values]}.

    posix=False so Windows paths keep their backslashes; the quotes it leaves
    attached are stripped by hand.
    """
    try:
        tokens = shlex.split(command, posix=False)
    except ValueError:
        return {}
    flags: dict[str, list[str]] = {}
    current = None
    for token in tokens[1:]:  # tokens[0] is the program name
        token = token.strip("\"'")
        if _is_flag(token):
            current = token
            flags.setdefault(current, [])
        elif current is not None:
            flags[current].append(token)
    return flags


def _flag_values(flags: dict[str, list[str]], *names: str) -> list[str]:
    """Values for the first of *names* present (flags have short and long forms)."""
    for name in names:
        if name in flags:
            return flags[name]
    return []


def describe_command(command: str) -> list[str]:
    """Summarize the parameters of a recorded run command."""
    flags = command_flags(command)
    lines = []

    queries = _flag_values(flags, "-q", "--queries")
    if queries:
        lines.append("Queries:   " + ", ".join(queries))

    dbs = _flag_values(flags, "-dbs", "--database")
    ns = _flag_values(flags, "-n", "--seqs")
    if dbs:
        labelled = [f"{db} (-n {ns[i]})" if i < len(ns) else db
                    for i, db in enumerate(dbs)]
        lines.append("Databases: " + ", ".join(labelled))

    opts = []
    for flag in ("--blast_type", "--aligner", "--tree_builder"):
        values = _flag_values(flags, flag)
        opts.append(values[0] if values else f"{CLI_DEFAULTS[flag]} (default)")
    lines.append("Options:   " + ", ".join(opts))

    extras = []
    for names, label in ((("-add", "--add_seqs"), "added seqs"),
                         (("-a", "--reroot"), "reroot"),
                         (("--motif", "--motifs"), "motifs"),
                         (("--hmm", "--hmms"), "HMMs"),
                         (("-aa", "--slice"), "slice"),
                         (("--datasets",), "datasets")):
        values = _flag_values(flags, *names)
        if values:
            extras.append(f"{label}: {' '.join(values)}")
    if extras:
        lines.append("Extras:    " + "; ".join(extras))
    return lines


def tree_stats(run_dir: Path) -> str:
    """Describe the combined tree: tip count, and how many tips per genome."""
    assets = run_dir / RUN_ASSETS_DIRNAME
    try:
        newick = (assets / COMBINED_TREE_NAME).read_text(encoding="utf-8",
                                                         errors="replace")
    except OSError:
        return ""
    # Tip labels follow '(' or ','; internal nodes carry support values, which
    # follow ')' instead and so are never captured.
    tips = [t.strip() for t in re.findall(r"[(,]\s*([^(),:;]+)", newick)]
    if not tips:
        return ""

    genomes: dict[str, str] = {}
    try:
        for line in (assets / MAPPING_NAME).read_text(
                encoding="utf-8", errors="replace").splitlines()[1:]:
            parts = line.split("\t")
            if len(parts) >= 2:
                genomes[parts[0]] = parts[1]
    except OSError:
        pass

    counts: dict[str, int] = {}
    for tip in tips:
        counts[genomes.get(tip, "unmapped")] = counts.get(genomes.get(tip, "unmapped"), 0) + 1
    breakdown = ", ".join(f"{g} {c}" for g, c in
                          sorted(counts.items(), key=lambda kv: -kv[1]))
    return f"{len(tips)} tips — {breakdown}"


def _data_rows(path: Path) -> list[list[str]]:
    """Rows of a pipeline TSV log, without its '#' preamble or header line."""
    try:
        lines = [ln for ln in path.read_text(encoding="utf-8",
                                             errors="replace").splitlines()
                 if ln.strip() and not ln.startswith("#")]
    except OSError:
        return []
    return [ln.split("\t") for ln in lines[1:]] if len(lines) > 1 else []


def log_stats(run_dir: Path) -> list[str]:
    """Summarize the de-duplication and translation logs, when present."""
    lines = []

    dedup = _data_rows(run_dir / DEDUP_LOG_NAME)
    if dedup:
        actions: dict[str, int] = {}
        for row in dedup:
            if len(row) > 4:
                actions[row[4]] = actions.get(row[4], 0) + 1
        summary = ", ".join(f"{c} {a}" for a, c in sorted(actions.items()))
        lines.append(f"De-dup:    {len(dedup)} entries — {summary}")

    report_path = run_dir / TRANSLATION_REPORT_NAME
    translated = _data_rows(report_path)
    if translated:
        stops = 0
        try:
            header = [ln for ln in report_path.read_text(
                encoding="utf-8", errors="replace").splitlines()
                if ln.strip() and not ln.startswith("#")][0].split("\t")
            col = header.index("n_internal_stops")
            stops = sum(1 for row in translated
                        if len(row) > col and row[col].strip() not in ("", "0"))
        except (OSError, ValueError, IndexError):
            pass
        lines.append(f"Translated: {len(translated)} sequences — "
                     f"{stops} with internal stops")
    return lines


def _is_wsl() -> bool:
    """Detect if running inside WSL."""
    try:
        return "microsoft" in Path("/proc/version").read_text().lower()
    except OSError:
        return False


def open_folder(path: Path):
    """Open a folder in the system file manager."""
    path = path.resolve()
    system = platform.system()
    if system == "Darwin":
        subprocess.Popen(["open", str(path)])
    elif system == "Windows":
        subprocess.Popen(["explorer", str(path)])
    else:
        if _is_wsl() and shutil.which("explorer.exe"):
            # Convert Linux path to Windows path for explorer.exe
            try:
                win_path = subprocess.check_output(
                    ["wslpath", "-w", str(path)], text=True).strip()
            except (subprocess.CalledProcessError, FileNotFoundError):
                win_path = str(path)
            subprocess.Popen(["explorer.exe", win_path])
        elif shutil.which("xdg-open"):
            subprocess.Popen(["xdg-open", str(path)])


def redraw_names(run_dir: Path, entry: str) -> list[str]:
    """Basenames of extra tree renders in a run, e.g. ENTRY_redraw_XI."""
    names = set()
    try:
        for f in run_dir.glob("*.pdf"):
            base = f.name[:-len(".pdf")]
            # Skip the derived companions (ENTRY.pdf.MSA.pdf and friends).
            if ".pdf" in base or base == entry:
                continue
            names.add(base)
    except OSError:
        pass
    return sorted(names)


class ToolTip:
    """Simple hover tooltip for a widget."""

    def __init__(self, widget, text):
        self.widget = widget
        self.text = text
        self.tip_window = None
        widget.bind("<Enter>", self._show)
        widget.bind("<Leave>", self._hide)

    def _show(self, event=None):
        if self.tip_window:
            return
        x = self.widget.winfo_rootx() + 20
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 2
        self.tip_window = tw = tk.Toplevel(self.widget)
        tw.wm_overrideredirect(True)
        tw.wm_geometry(f"+{x}+{y}")
        label = tk.Label(tw, text=self.text, background="#ffffe0",
                         relief="solid", borderwidth=1, font=("TkDefaultFont", 9))
        label.pack()

    def _hide(self, event=None):
        if self.tip_window:
            self.tip_window.destroy()
            self.tip_window = None


class GenomeSelectorApp:
    def __init__(self, root: tk.Tk, default_n: int | None = None,
                 selected_genomes: list[str] | None = None,
                 restore: bool = True):
        self.root = root
        root.title("BAT Genome Selector")
        # Tall enough that the genome table and the Recent Runs details
        # panel are both usable; below this one of them collapses.
        root.minsize(1050, 800)

        self.default_n = DEFAULT_N if default_n is None else default_n
        self._default_n_from_cli = default_n is not None
        self.selected_genomes = {
            name.replace("\\", "/").casefold()
            for name in (selected_genomes or [])
        }
        self.row_tokens: dict[str, list[str]] = {}

        self._build_banner()
        # Before the genome table on purpose: pack hands out height in the order
        # widgets are packed, and clips whatever is packed last. The genome list
        # scrolls, so it can absorb a short window; the run details panel cannot,
        # and used to be squeezed to nothing.
        self._build_output_panel()
        self._build_genome_table()
        self._build_controls()
        self._build_options_panel()
        self._build_advanced_panel()
        self._build_action_buttons()

        # Everything above is the state this launcher starts in — command-line
        # options included. "Reset All" restores exactly this.
        self._launch_defaults = self._snapshot()
        if restore:
            self._load_settings()
        root.protocol("WM_DELETE_WINDOW", self._on_close)

        # Fix Ctrl+A select-all in Entry/Combobox/Spinbox widgets
        def _select_all_text(event):
            event.widget.select_range(0, "end")
            event.widget.icursor("end")
            return "break"

        for cls in (ttk.Entry, ttk.Combobox, ttk.Spinbox):
            root.bind_class(cls, "<Control-a>", _select_all_text)

        # Sync header widths once all widgets have been laid out
        root.after_idle(lambda: (root.update_idletasks(), self._sync_columns()))

    # ------------------------------------------------------------------
    # Banner
    # ------------------------------------------------------------------
    def _build_banner(self):
        bat_col = 45
        rows = [
            ("        ,-- ATGCATGC  H. sapiens",       ""),
            ("    ,---+",                              "     _.-~^~-._        /\\  /\\        _.-~^~-._"),
            ("    |   `-- ATGCGTGC  D. melanogaster",  "   ,'         `-.____/    \\____.-'`         `,"),
            (" ---+",                                  "  /                 B.A.T.                  \\"),
            ("    |   ,-- ATGCATCC  A. thaliana",      "  \\__,           _,-~      ~-,_          ,__/"),
            ("    `---+",                              "      `-.____.-'                `-.____.-'"),
            ("        `-- ATGCGTCC  O. sativa",        ""),
        ]
        body = "\n".join(
            (tree + " " * max(1, bat_col - len(tree)) + bat) if bat else tree
            for tree, bat in rows
        )
        banner = (
            "  B.A.T.  BLAST - ALIGN - TREE\n"
            "  -------------------------------------------------------------------------------------\n"
            + body
        )
        tk.Label(self.root, text=banner, font=("Courier", 9), justify="left",
                 anchor="w", padx=12, pady=4).pack(fill="x")

        exts = "  ".join(sorted(FASTA_EXTENSIONS))
        summary = (
            f"Scanning genomes/ and subfolders for FASTA files ({exts})\n"
            "1. Fill in query IDs, select hit databases to search, and set -hdr tokens to parse descriptions to the parts you want. \n"
            "2. Click 'Generate Command', then 'Copy to Clipboard'.\n"
            "3. Paste and run the command in a terminal with your conda/mamba environment activated."
        )
        tk.Label(self.root, text=summary, font=("TkDefaultFont", 9),
                 justify="left", anchor="w", padx=12).pack(fill="x")

    # ------------------------------------------------------------------
    # Genome table (scrollable)
    # ------------------------------------------------------------------
    def _build_genome_table(self):
        # Fixed column headers
        self.header_frame = header_frame = tk.Frame(self.root)
        header_frame.pack(fill="x", padx=8, pady=(8, 0))

        # Section headings (row 0)
        # "Query" spans cols 0-2, separator at col 3, "Database Settings" spans cols 4-10
        query_heading = tk.Label(header_frame, text="Query",
                                 font=("TkDefaultFont", 9, "bold", "italic"),
                                 fg="#444")
        query_heading.grid(row=0, column=0, columnspan=3, sticky="ew", padx=4)
        db_heading = tk.Label(header_frame, text="Database Settings",
                              font=("TkDefaultFont", 9, "bold", "italic"),
                              fg="#444")
        db_heading.grid(row=0, column=4, columnspan=7, sticky="ew", padx=4)

        # Column headers (row 1) - col 3 reserved for vertical divider
        for col, (text, px) in [
            (0, ("Genome file", 4)),   (1, ("Ex.", 2)),
            (2, ("-q (queries)", 4)),
            # col 3 = divider
            (4, ("Include", (4, 2))),  (5, ("-hdr", 4)),
            (6, ("Parsed name", 4)),
            (7, ("-hdr_sfx", 4)),
            (8, ("-n", 4)),            (9, ("Type", 4)),
            (10, ("DB", 4)),
        ]:
            lbl = tk.Label(header_frame, text=text, font=("TkDefaultFont", 9, "bold"),
                           anchor="center")
            lbl.grid(row=1, column=col, padx=px, sticky="ew")

        # Vertical divider between query and database sections
        header_frame.columnconfigure(DIVIDER_COL, minsize=DIVIDER_COL_WIDTH, weight=0)
        ttk.Separator(header_frame, orient="vertical").grid(
            row=0, column=DIVIDER_COL, rowspan=2, sticky="ns")

        # Scrollable genome list
        list_frame = tk.Frame(self.root)
        list_frame.pack(fill="both", expand=True, padx=8)

        self.canvas = canvas = tk.Canvas(list_frame, highlightthickness=0)
        scrollbar = ttk.Scrollbar(list_frame, orient="vertical", command=canvas.yview)
        self.inner_frame = tk.Frame(canvas)
        self.inner_frame.columnconfigure(DIVIDER_COL, minsize=DIVIDER_COL_WIDTH, weight=0)

        self.inner_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        canvas.create_window((0, 0), window=self.inner_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        # Mousewheel scrolling
        def _on_mousewheel(event):
            canvas.yview_scroll(-1 * (event.delta // 120 or (
                -1 if event.num == 5 else 1)), "units")

        canvas.bind_all("<MouseWheel>", _on_mousewheel)
        canvas.bind_all("<Button-4>", _on_mousewheel)
        canvas.bind_all("<Button-5>", _on_mousewheel)

        # Update scroll region whenever the inner frame resizes. Column
        # syncing is handled explicitly on init / refresh to avoid recursive
        # Configure events (sync resizes columns → triggers Configure → loop).
        self.inner_frame.bind("<Configure>",
                              lambda e: canvas.configure(scrollregion=canvas.bbox("all")))

        self._populate_rows()

    # ------------------------------------------------------------------
    # Control buttons (select/deselect/clear/refresh)
    # ------------------------------------------------------------------
    def _build_controls(self):
        btn_frame = tk.Frame(self.root)
        btn_frame.pack(fill="x", padx=8, pady=4)

        ttk.Button(btn_frame, text="Select All Hits",
                   command=self._select_all).pack(side="left", padx=(0, 4))
        ttk.Button(btn_frame, text="Deselect All Hits",
                   command=self._deselect_all).pack(side="left", padx=(0, 4))
        ttk.Button(btn_frame, text="Clear Fields",
                   command=self._clear_fields).pack(side="left", padx=(0, 4))
        ttk.Button(btn_frame, text="Refresh",
                   command=self._refresh).pack(side="left", padx=(0, 8))

        ttk.Separator(btn_frame, orient="vertical").pack(side="left", fill="y", padx=(0, 8))
        tk.Label(btn_frame, text="Default -n:").pack(side="left")
        self.default_n_var = tk.StringVar(value=str(self.default_n))
        default_n_spin = ttk.Spinbox(btn_frame, textvariable=self.default_n_var,
                                     from_=MIN_N, to=MAX_N, width=5,
                                     command=self._set_default_n)
        default_n_spin.pack(side="left", padx=(2, 4))
        default_n_spin.bind("<Return>", lambda _e: self._set_default_n())
        ToolTip(default_n_spin, "Default number of hits for new rows and Clear Fields")
        ttk.Button(btn_frame, text="Set All -n",
                   command=self._set_all_n).pack(side="left")

        # Kept away from the buttons used routinely — this one discards work.
        reset_btn = ttk.Button(btn_frame, text="Reset All",
                               command=self._reset_all)
        reset_btn.pack(side="right")
        ToolTip(reset_btn,
                "Clear every field back to the state this launcher starts in")
        ttk.Separator(btn_frame, orient="vertical").pack(
            side="right", fill="y", padx=8)

        # Selected hit databases display
        self.selected_label = tk.Label(self.root, text="Hit DBs: (none)",
                                       font=("TkDefaultFont", 9), anchor="w",
                                       justify="left", wraplength=1000)
        self.selected_label.pack(fill="x", padx=12, pady=(0, 4))

    # ------------------------------------------------------------------
    # Command options (aligner, tree, blast_type, threads)
    # ------------------------------------------------------------------
    def _build_options_panel(self):
        opt_frame = tk.Frame(self.root)
        opt_frame.pack(fill="x", padx=8, pady=(0, 4))

        tk.Label(opt_frame, text="Aligner:").pack(side="left")
        self.aligner_var = tk.StringVar(value="clustalo")
        aligner_combo = ttk.Combobox(opt_frame, textvariable=self.aligner_var,
                                     values=["clustalo", "mafft"], width=9,
                                     state="readonly")
        aligner_combo.pack(side="left", padx=(2, 8))

        tk.Label(opt_frame, text="MAFFT mode:").pack(side="left")
        self.mafft_mode_var = tk.StringVar(value="auto")
        self.mafft_combo = ttk.Combobox(opt_frame, textvariable=self.mafft_mode_var,
                                        values=["auto", "linsi", "einsi", "fftns2"],
                                        width=7, state="disabled")
        self.mafft_combo.pack(side="left", padx=(2, 8))

        def _on_aligner_change(*_args):
            if self.aligner_var.get() == "mafft":
                self.mafft_combo.configure(state="readonly")
            else:
                self.mafft_combo.configure(state="disabled")

        self.aligner_var.trace_add("write", _on_aligner_change)

        tk.Label(opt_frame, text="Tree builder:").pack(side="left")
        self.tree_var = tk.StringVar(value="FastTree")
        ttk.Combobox(opt_frame, textvariable=self.tree_var,
                     values=["FastTree", "RAxML"], width=9,
                     state="readonly").pack(side="left", padx=(2, 8))

        tk.Label(opt_frame, text="BLAST type:").pack(side="left")
        self.blast_type_var = tk.StringVar(value="tblastn")
        ttk.Combobox(opt_frame, textvariable=self.blast_type_var,
                     values=["tblastn", "blastp"], width=8,
                     state="readonly").pack(side="left", padx=(2, 8))

        tk.Label(opt_frame, text="Threads:").pack(side="left")
        self.threads_var = tk.StringVar(value=str(max(1, os.cpu_count() // 2)))
        ttk.Spinbox(opt_frame, textvariable=self.threads_var,
                    from_=1, to=os.cpu_count() or 8, width=4).pack(side="left", padx=(2, 0))

    # ------------------------------------------------------------------
    # Advanced options (collapsible): outgroups, AA slice, motifs, HMM
    # ------------------------------------------------------------------
    def _build_advanced_panel(self):
        self._adv_visible = tk.BooleanVar(value=False)
        toggle_frame = tk.Frame(self.root)
        toggle_frame.pack(fill="x", padx=8, pady=(0, 2))
        self._adv_toggle_btn = ttk.Button(
            toggle_frame, text="> Advanced Options",
            command=self._toggle_advanced)
        self._adv_toggle_btn.pack(side="left")

        self.adv_frame = tk.Frame(self.root)
        # Not packed yet — toggled on demand

        # Row 1: Outgroups
        row1 = tk.Frame(self.adv_frame)
        row1.pack(fill="x", pady=2)
        tk.Label(row1, text="-add (outgroup seqs):").pack(side="left")
        self.add_seqs_var = tk.StringVar()
        ttk.Entry(row1, textvariable=self.add_seqs_var, width=30).pack(side="left", padx=(2, 8))
        ToolTip(row1.winfo_children()[-1], "Space-separated sequence IDs to add as outgroups")

        tk.Label(row1, text="-add_db (outgroup DBs):").pack(side="left")
        self.add_dbs_var = tk.StringVar()
        ttk.Entry(row1, textvariable=self.add_dbs_var, width=30).pack(side="left", padx=(2, 8))
        ToolTip(row1.winfo_children()[-1], "Space-separated FASTA files for outgroup sequences")

        tk.Label(row1, text="-a (reroot on):").pack(side="left")
        self.reroot_var = tk.StringVar()
        ttk.Entry(row1, textvariable=self.reroot_var, width=18).pack(side="left", padx=(2, 0))
        ToolTip(row1.winfo_children()[-1],
                "One tip to root the tree on, usually an outgroup added above.\n"
                "Use the ID as it appears in the tree, i.e. after -hdr parsing\n"
                "(e.g. AT2G38240, not AT2G38240.1)")

        # Row 2: AA slice
        row2 = tk.Frame(self.adv_frame)
        row2.pack(fill="x", pady=2)
        tk.Label(row2, text="-aa / --slice (AA range):").pack(side="left")
        self.aa_start_var = tk.StringVar()
        ttk.Entry(row2, textvariable=self.aa_start_var, width=6).pack(side="left", padx=(2, 2))
        tk.Label(row2, text="to").pack(side="left")
        self.aa_end_var = tk.StringVar()
        ttk.Entry(row2, textvariable=self.aa_end_var, width=6).pack(side="left", padx=(2, 8))
        ToolTip(row2.winfo_children()[1], "Start AA position (e.g. 10)")
        ToolTip(row2.winfo_children()[3], "End AA position (e.g. 200)")

        # Row 3: Motifs
        row3 = tk.Frame(self.adv_frame)
        row3.pack(fill="x", pady=2)
        tk.Label(row3, text="--motif:").pack(side="left")
        self.motif_var = tk.StringVar()
        ttk.Entry(row3, textvariable=self.motif_var, width=40).pack(side="left", padx=(2, 8))
        ToolTip(row3.winfo_children()[-1],
                "Space-separated motif patterns (NAME=PATTERN or just PATTERN)")

        tk.Label(row3, text="Syntax:").pack(side="left")
        self.motif_syntax_var = tk.StringVar(value="regex")
        ttk.Combobox(row3, textvariable=self.motif_syntax_var,
                     values=["regex", "prosite"], width=8,
                     state="readonly").pack(side="left", padx=(2, 8))

        self.motif_overlap_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(row3, text="Allow overlap",
                        variable=self.motif_overlap_var).pack(side="left")

        # Row 4: HMM
        row4 = tk.Frame(self.adv_frame)
        row4.pack(fill="x", pady=2)
        tk.Label(row4, text="--hmm:").pack(side="left")
        self.hmm_var = tk.StringVar()
        ttk.Entry(row4, textvariable=self.hmm_var, width=40).pack(side="left", padx=(2, 0))
        ToolTip(row4.winfo_children()[-1], "Space-separated HMMER profile files (.hmm)")

    def _toggle_advanced(self):
        if self._adv_visible.get():
            self.adv_frame.pack_forget()
            self._adv_toggle_btn.configure(text="> Advanced Options")
            self._adv_visible.set(False)
        else:
            # Insert advanced frame just after the toggle button's parent
            self.adv_frame.pack(fill="x", padx=20, pady=(0, 4),
                                after=self._adv_toggle_btn.master)
            self._adv_toggle_btn.configure(text="v Advanced Options")
            self._adv_visible.set(True)

    # ------------------------------------------------------------------
    # Action buttons (Generate, Copy)
    # ------------------------------------------------------------------
    def _build_action_buttons(self):
        gen_frame = tk.Frame(self.root)
        gen_frame.pack(fill="x", padx=8, pady=(0, 4))

        self.prepend_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(gen_frame, text="Prepend  blast-align-tree",
                        variable=self.prepend_var).pack(side="left", padx=(0, 8))
        tk.Button(gen_frame, text="Generate Command",
                  command=self._generate,
                  font=("TkDefaultFont", 9, "bold"),
                  bg="#ffd166", activebackground="#f4b942",
                  fg="#111111", activeforeground="#111111",
                  relief="raised", borderwidth=2).pack(side="left", padx=(0, 4))
        ttk.Button(gen_frame, text="Copy to Clipboard",
                   command=self._copy_to_clipboard).pack(side="left", padx=(0, 4))

    # ------------------------------------------------------------------
    # Output panel (tabbed: Command | Recent Runs)
    # ------------------------------------------------------------------
    def _build_output_panel(self):
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(side="bottom", fill="both", expand=True,
                           padx=8, pady=(0, 8))

        # Command tab
        cmd_frame = tk.Frame(self.notebook)
        self.output = tk.Text(cmd_frame, height=6, wrap="none",
                              font=("Courier", 10))
        self.output.bind("<Key>", lambda e: "break" if e.keysym not in
                         ("c", "a") or not (e.state & 0x4) else None)
        self.output.pack(fill="both", expand=True)
        self.notebook.add(cmd_frame, text="Command")

        # Recent Runs tab
        self._build_recent_runs_tab()

    def _build_recent_runs_tab(self):
        runs_frame = tk.Frame(self.notebook)

        # Toolbar
        toolbar = tk.Frame(runs_frame)
        toolbar.pack(fill="x", pady=(4, 2), padx=4)
        ttk.Button(toolbar, text="Refresh", command=self._refresh_runs).pack(side="left")

        # Packed before the list: the list canvas expands, so anything packed
        # after it is squeezed out of the tab entirely.
        details = self._build_details_panel(runs_frame)
        details.pack(side="bottom", fill="x", padx=4, pady=(0, 6))

        # Scrollable list
        list_frame = tk.Frame(runs_frame)
        list_frame.pack(side="top", fill="both", expand=True, padx=4, pady=(0, 4))

        self.runs_canvas = tk.Canvas(list_frame, highlightthickness=0,
                                     height=RUNS_LIST_HEIGHT)
        runs_scroll = ttk.Scrollbar(list_frame, orient="vertical",
                                    command=self.runs_canvas.yview)
        self.runs_inner = tk.Frame(self.runs_canvas)
        self.runs_inner.bind(
            "<Configure>",
            lambda e: self.runs_canvas.configure(scrollregion=self.runs_canvas.bbox("all")))
        self.runs_canvas.create_window((0, 0), window=self.runs_inner, anchor="nw")
        self.runs_canvas.configure(yscrollcommand=runs_scroll.set)
        self.runs_canvas.pack(side="left", fill="both", expand=True)
        runs_scroll.pack(side="right", fill="y")

        self.notebook.add(runs_frame, text="Recent Runs")
        self._refresh_runs()

    def _build_details_panel(self, parent) -> tk.LabelFrame:
        """Details for the selected run: what it was built from, what came out,
        the command that produced it, and the hint for re-drawing its trees
        (the same one the pipeline prints when it finishes). Defaults to the
        newest run; any row's buttons swap in that run and copy from it.

        One read-only text rather than a stack of boxes: the notebook only gets
        the height the panels above it leave over, and every line spent here is
        a run the list cannot show.
        """
        details = tk.LabelFrame(parent, text="Run details",
                                font=("TkDefaultFont", 9, "bold"))

        self.runs_info = tk.Text(details, height=DETAILS_MAX_LINES, wrap="none",
                                 font=("Courier", 9), relief="flat",
                                 background=self.root.cget("background"))
        self.runs_info.bind("<Key>", lambda e: "break" if e.keysym not in
                            ("c", "a") or not (e.state & 0x4) else None)
        # Commands run past the window; scroll rather than wrap, so the aligned
        # detail lines above them stay aligned.
        info_scroll = ttk.Scrollbar(details, orient="horizontal",
                                    command=self.runs_info.xview)
        self.runs_info.configure(xscrollcommand=info_scroll.set)
        self.runs_info.pack(fill="x", padx=4, pady=(2, 0))
        info_scroll.pack(fill="x", padx=4, pady=(0, 4))
        return details

    def _refresh_runs(self):
        """Populate the Recent Runs list."""
        for widget in self.runs_inner.winfo_children():
            widget.destroy()

        runs = scan_recent_runs()
        if not runs:
            tk.Label(self.runs_inner, text="No runs found.",
                     font=("TkDefaultFont", 9), fg="gray").grid(
                row=0, column=0, padx=8, pady=8)
            self._show_run(None)
            return

        # Column headers
        for col, text in enumerate(["Entry", "Timestamp", "Contents", "", "", "", ""]):
            tk.Label(self.runs_inner, text=text,
                     font=("TkDefaultFont", 9, "bold")).grid(
                row=0, column=col, padx=6, pady=(4, 2), sticky="w")

        for i, (entry, timestamp, path) in enumerate(runs, start=1):
            display_ts = format_timestamp(timestamp)

            tk.Label(self.runs_inner, text=entry, anchor="w").grid(
                row=i, column=0, padx=6, pady=1, sticky="w")
            tk.Label(self.runs_inner, text=display_ts, anchor="w").grid(
                row=i, column=1, padx=6, pady=1, sticky="w")

            # Summarize report files from the run folder
            contents = self._summarize_run(path)
            tk.Label(self.runs_inner, text=contents, anchor="w",
                     font=("TkDefaultFont", 8), fg="#555").grid(
                row=i, column=2, padx=6, pady=1, sticky="w")

            ttk.Button(self.runs_inner, text="Open Run", width=11,
                       command=lambda p=path: open_folder(p)).grid(
                row=i, column=3, padx=6, pady=1)

            run = (entry, timestamp, path)
            ttk.Button(self.runs_inner, text="Details", width=9,
                       command=lambda r=run: self._show_run(r)).grid(
                row=i, column=4, padx=(0, 4), pady=1)

            # Copying selects the run too, so the panel always shows what was
            # just put on the clipboard.
            rerun = ttk.Button(
                self.runs_inner, text="Re-run", width=9,
                command=lambda r=run: self._copy_for_run(r, "run_command"))
            rerun.grid(row=i, column=5, padx=(0, 4), pady=1)
            if read_run_command(path) is None:
                # Nothing to copy: this run predates command logging.
                rerun.state(["disabled"])
                ToolTip(rerun, "This run was archived before blast-align-tree "
                               "started saving its command")

            redraw = ttk.Button(
                self.runs_inner, text="Re-draw", width=9,
                command=lambda r=run: self._copy_for_run(r, "redraw_cmd"))
            redraw.grid(row=i, column=6, padx=(0, 6), pady=1)
            ToolTip(redraw, "Copy the Rscript command that re-draws this run's "
                            "trees, e.g. from a subnode")

        self._show_run(runs[0])

    def _copy_for_run(self, run: tuple[str, str, Path], which: str):
        """Select a run, then copy one of its two commands."""
        self._show_run(run)
        self._copy(getattr(self, which))

    def _show_run(self, run: tuple[str, str, Path] | None):
        """Fill the details panel from a (entry, timestamp, path) run."""
        if run is None:
            self.run_command = None
            self.redraw_cmd = None
            self._set_text(self.runs_info, "")
            return

        entry, timestamp, path = run
        command = read_run_command(path)
        self.run_command = command
        datasets = command_flags(command).get("--datasets") if command else None
        self.redraw_cmd = redraw_command(entry, timestamp,
                                         datasets[0] if datasets else None)

        lines = [f"Run:       {entry}   {format_timestamp(timestamp)}"]
        if command:
            lines.extend(describe_command(command))
        else:
            lines.append("Params:    not recorded (run predates command logging)")
        tree = tree_stats(path)
        if tree:
            lines.append(f"Tree:      {tree}")
        lines.extend(log_stats(path))
        outputs = self._summarize_run(path)
        if outputs:
            lines.append(f"Outputs:   {outputs}")
        redraws = redraw_names(path, entry)
        if redraws:
            lines.append(f"Re-draws:  {', '.join(redraws)}")

        lines.append("")
        lines.append("Re-run:    " + (command or
                     "not recorded — archived before blast-align-tree "
                     "started saving its command"))
        lines.append(f"Re-draw:   {self.redraw_cmd}")
        self._set_text(self.runs_info, "\n".join(lines))
        # Only as tall as it needs to be; the run list gets what is left.
        self.runs_info.configure(height=min(len(lines), DETAILS_MAX_LINES))

    @staticmethod
    def _set_text(widget: tk.Text, text: str):
        widget.delete("1.0", "end")
        widget.insert("1.0", text)

    def _copy(self, text: str | None):
        if text:
            self.root.clipboard_clear()
            self.root.clipboard_append(text)
            self.root.update()


    @staticmethod
    def _summarize_run(path: Path) -> str:
        """Return a short summary of file types in a run directory."""
        counts: dict[str, int] = {}
        try:
            for f in path.iterdir():
                if f.is_file():
                    ext = f.suffix.lower()
                    if ext in (".pdf",):
                        counts["PDFs"] = counts.get("PDFs", 0) + 1
                    elif ext in (".nwk", ".newick", ".tree"):
                        counts["trees"] = counts.get("trees", 0) + 1
                    elif ext in (".fa", ".fasta", ".faa"):
                        counts["FASTAs"] = counts.get("FASTAs", 0) + 1
        except OSError:
            pass
        if not counts:
            return ""
        return ", ".join(f"{v} {k}" for k, v in counts.items())

    # ------------------------------------------------------------------
    # Column sync
    # ------------------------------------------------------------------
    def _sync_columns(self):
        """Align header_frame and inner_frame columns by forcing each matching
        column in both frames to the larger of the two natural cell widths."""
        num_cols = max(self.inner_frame.grid_size()[0],
                       self.header_frame.grid_size()[0])
        # Reset minsize so grid_bbox reports the true natural widths.
        for col in range(num_cols):
            if col == DIVIDER_COL:
                continue
            self.inner_frame.columnconfigure(col, minsize=0)
            self.header_frame.columnconfigure(col, minsize=0)
        self.root.update_idletasks()
        # grid_bbox returns (x, y, w, h) with w including padx — use it so
        # padding contributes to alignment, not just widget width.
        for col in range(num_cols):
            if col == DIVIDER_COL:
                continue
            try:
                inner_bbox = self.inner_frame.grid_bbox(column=col, row=0)
                header_bbox = self.header_frame.grid_bbox(column=col, row=1)
            except tk.TclError:
                continue
            inner_w = inner_bbox[2] if inner_bbox else 0
            header_w = header_bbox[2] if header_bbox else 0
            min_w = max(inner_w, header_w)
            if min_w:
                self.inner_frame.columnconfigure(col, minsize=min_w)
                self.header_frame.columnconfigure(col, minsize=min_w)
        # Keep the divider column fixed in both frames so the vertical rules align.
        self.header_frame.columnconfigure(DIVIDER_COL, minsize=DIVIDER_COL_WIDTH, weight=0)
        self.inner_frame.columnconfigure(DIVIDER_COL, minsize=DIVIDER_COL_WIDTH, weight=0)

    # ------------------------------------------------------------------
    # Genome rows
    # ------------------------------------------------------------------
    def _populate_rows(self):
        """Build genome rows in the scrollable inner_frame."""
        genome_files = scan_genomes()
        self.rows = []
        self.row_tokens = {}

        for i, gpath in enumerate(genome_files):
            tokens = detect_header_tokens(gpath)
            if "id" in tokens:
                tokens.remove("id")
            tokens.insert(0, "id")
            default_hdr = "id"

            rel = gpath.relative_to(GENOMES_DIR).as_posix()
            self.row_tokens[rel] = tokens
            chk_var = tk.BooleanVar(value=rel.casefold() in self.selected_genomes)
            hdr_var = tk.StringVar(value=default_hdr)
            sfx_var = tk.StringVar(value="none")
            n_var = tk.StringVar(value=str(self.default_n))
            q_var = tk.StringVar(value="")
            example_header = get_example_header(gpath)
            preview_var = tk.StringVar(
                value=parse_header_token(example_header, default_hdr, sfx_var.get()))

            display_name = truncate_name(rel)
            lbl = tk.Label(self.inner_frame, text=display_name, anchor="w")
            lbl.grid(row=i, column=0, sticky="w", padx=4)
            if len(rel) > MAX_LABEL_LEN:
                ToolTip(lbl, rel)

            ttk.Button(
                self.inner_frame, text="Ex.", width=3,
                command=lambda qv=q_var, gp=gpath: self._fill_example(qv, gp)
            ).grid(row=i, column=1, padx=2, pady=1)

            ttk.Entry(self.inner_frame, textvariable=q_var, width=20).grid(
                row=i, column=2, sticky="w", padx=4, pady=1)

            # Vertical divider between query and database sections
            ttk.Separator(self.inner_frame, orient="vertical").grid(
                row=i, column=DIVIDER_COL, sticky="ns")

            ttk.Checkbutton(self.inner_frame, variable=chk_var).grid(
                row=i, column=4, padx=(4, 2))

            ttk.Combobox(self.inner_frame, textvariable=hdr_var,
                         values=tokens, width=18).grid(
                row=i, column=5, padx=4, pady=1)

            preview_entry = ttk.Entry(self.inner_frame, textvariable=preview_var,
                                      width=22, state="readonly")
            preview_entry.grid(row=i, column=6, padx=4, pady=1)

            def _update_preview(*_args, hv=hdr_var, sv=sfx_var, pv=preview_var,
                                eh=example_header):
                pv.set(parse_header_token(eh, hv.get(), sv.get()))

            hdr_var.trace_add("write", _update_preview)
            sfx_var.trace_add("write", _update_preview)

            ttk.Entry(self.inner_frame, textvariable=sfx_var, width=8).grid(
                row=i, column=7, padx=4, pady=1)

            ttk.Spinbox(self.inner_frame, textvariable=n_var,
                        from_=MIN_N, to=MAX_N, width=5).grid(
                row=i, column=8, padx=4, pady=1)

            tk.Label(self.inner_frame, text=dbtype_for(gpath)).grid(
                row=i, column=9, padx=4, pady=1)

            if has_blastdb(gpath):
                tk.Label(self.inner_frame, text="Yes", fg="green").grid(
                    row=i, column=10, padx=4, pady=1)
            else:
                ttk.Button(
                    self.inner_frame, text="makedb", width=6,
                    command=lambda gp=gpath: self._copy_makeblastdb_single(gp)
                ).grid(row=i, column=10, padx=4, pady=1)

            chk_var.trace_add("write", self._update_selected_display)
            self.rows.append((chk_var, hdr_var, sfx_var, n_var, q_var, gpath, rel))

    # ------------------------------------------------------------------
    # Selection helpers
    # ------------------------------------------------------------------
    def _update_selected_display(self, *_args):
        names = [name for chk_var, _h, _s, _n, _q, _g, name in self.rows
                 if chk_var.get()]
        if names:
            self.selected_label.configure(
                text=f"Hit DBs ({len(names)}): " + ", ".join(names))
        else:
            self.selected_label.configure(text="Hit DBs: (none)")

    # ------------------------------------------------------------------
    # Saved settings (snapshot / restore / reset)
    # ------------------------------------------------------------------
    def _snapshot_rows(self) -> dict:
        return {
            rel: {
                "selected": chk_var.get(),
                "hdr": hdr_var.get(),
                "sfx": sfx_var.get(),
                "n": n_var.get(),
                "q": q_var.get(),
            }
            for chk_var, hdr_var, sfx_var, n_var, q_var, _gpath, rel in self.rows
        }

    def _snapshot(self) -> dict:
        globals_ = {key: getattr(self, attr).get()
                    for key, attr in GLOBAL_SETTING_VARS.items()}
        globals_["adv_visible"] = self._adv_visible.get()
        return {
            "version": SETTINGS_VERSION,
            "globals": globals_,
            "rows": self._snapshot_rows(),
        }

    def _apply_rows(self, saved_rows: dict) -> list[str]:
        """Apply saved per-row state. Returns warnings for stale headers."""
        warnings = []
        for chk_var, hdr_var, sfx_var, n_var, q_var, _gpath, rel in self.rows:
            row = saved_rows.get(rel)
            if not isinstance(row, dict):
                continue  # genome added since the snapshot — keep its defaults
            chk_var.set(bool(row.get("selected", False)))

            hdr = str(row.get("hdr", "id"))
            tokens = self.row_tokens.get(rel, [])
            if tokens and hdr not in tokens:
                warnings.append(
                    f"  {rel}: header '{hdr}' is no longer in this file — using 'id'")
                hdr = "id"
            hdr_var.set(hdr)

            sfx_var.set(str(row.get("sfx", "none")))
            n_var.set(str(row.get("n", self.default_n)))
            q_var.set(str(row.get("q", "")))
        return warnings

    def _apply_settings(self, data: dict) -> list[str]:
        globals_ = data.get("globals") or {}
        for key, attr in GLOBAL_SETTING_VARS.items():
            if key not in globals_:
                continue
            var = getattr(self, attr)
            value = globals_[key]
            var.set(bool(value) if isinstance(var, tk.BooleanVar) else str(value))
        # Keep self.default_n (used for new rows) in step with the spinbox.
        self._set_default_n()

        if bool(globals_.get("adv_visible", False)) != self._adv_visible.get():
            self._toggle_advanced()

        return self._apply_rows(data.get("rows") or {})

    def _save_settings(self):
        try:
            SETTINGS_FILE.write_text(
                json.dumps(self._snapshot(), indent=2), encoding="utf-8")
        except OSError as exc:
            self._set_output(f"Could not save settings to {SETTINGS_FILE.name}: {exc}")

    def _load_settings(self):
        if not SETTINGS_FILE.exists():
            return
        try:
            data = json.loads(SETTINGS_FILE.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            self._set_output(f"Ignoring unreadable {SETTINGS_FILE.name}: {exc}")
            return
        if not isinstance(data, dict) or data.get("version") != SETTINGS_VERSION:
            self._set_output(
                f"Ignoring {SETTINGS_FILE.name}: unrecognized settings format.")
            return

        # Command-line options win over the saved session.
        if self._default_n_from_cli:
            (data.get("globals") or {}).pop("default_n", None)

        warnings = self._apply_settings(data)

        # Genomes named on the command line win over the saved selection.
        if self.selected_genomes:
            for chk_var, _h, _s, _n, _q, _gpath, rel in self.rows:
                chk_var.set(rel.casefold() in self.selected_genomes)

        if warnings:
            self._set_output("Restored previous settings, with changes:\n"
                             + "\n".join(warnings))

    def _reset_all(self):
        confirm = messagebox.askyesno(
            "Reset All",
            "Reset every field to the state this launcher starts in?\n\n"
            "Queries, header choices and advanced options will be cleared.",
            default=messagebox.NO,
            parent=self.root,
        )
        if not confirm:
            return
        self._apply_settings(self._launch_defaults)
        self._save_settings()
        self._set_output("")

    def _on_close(self):
        self._save_settings()
        self.root.destroy()

    def _set_default_n(self) -> int | None:
        raw = self.default_n_var.get().strip()
        try:
            self.default_n = positive_int(raw)
        except argparse.ArgumentTypeError as exc:
            self.default_n_var.set(str(self.default_n))
            self._set_output(f"Default -n {exc}.")
            return None
        self.default_n_var.set(str(self.default_n))
        return self.default_n

    def _set_all_n(self):
        default_n = self._set_default_n()
        if default_n is None:
            return
        for _chk_var, _hdr_var, _sfx_var, n_var, _q_var, _gpath, _name in self.rows:
            n_var.set(str(default_n))

    def _refresh(self):
        if hasattr(self, "default_n_var"):
            self._set_default_n()
        # Rescanning rebuilds every row from scratch, so carry the current
        # per-row state across; genomes added since are left at their defaults.
        previous = self._snapshot_rows() if getattr(self, "rows", None) else {}
        _id_cache.clear()
        for widget in self.inner_frame.winfo_children():
            widget.destroy()
        self._populate_rows()
        warnings = self._apply_rows(previous)
        self.root.update_idletasks()
        self._sync_columns()
        self._update_selected_display()
        if warnings:
            self._set_output("Refreshed, with changes:\n" + "\n".join(warnings))

    def _select_all(self):
        for chk_var, *_ in self.rows:
            chk_var.set(True)

    def _deselect_all(self):
        for chk_var, *_ in self.rows:
            chk_var.set(False)

    def _fill_example(self, q_var: tk.StringVar, gpath: Path):
        ids = load_fasta_ids(gpath)
        if ids:
            q_var.set(random.choice(list(ids)))

    def _copy_makeblastdb_single(self, gpath: Path):
        dt = dbtype_for(gpath)
        rel = gpath.relative_to(GENOMES_DIR).as_posix()
        cmd = f"makeblastdb -in genomes/{rel} -dbtype {dt} -parse_seqids"
        self._set_output(cmd)
        self.root.clipboard_clear()
        self.root.clipboard_append(cmd)

    def _clear_fields(self):
        default_n = self._set_default_n()
        if default_n is None:
            return
        for chk_var, hdr_var, sfx_var, n_var, q_var, gpath, name in self.rows:
            q_var.set("")
            sfx_var.set("none")
            n_var.set(str(default_n))

    # ------------------------------------------------------------------
    # Command generation
    # ------------------------------------------------------------------
    @staticmethod
    def _parse_queries(raw: str) -> list[str]:
        return [q for q in re.split(r'[,\s]+', raw.strip()) if q]

    def _build_args(self) -> list[str] | None:
        """Build the argument list for blast-align-tree. Returns None on error."""
        hit_dbs = []
        query_rows = []
        for chk_var, hdr_var, sfx_var, n_var, q_var, gpath, name in self.rows:
            if chk_var.get():
                hit_dbs.append((name, hdr_var.get(), sfx_var.get().strip(),
                                n_var.get(), gpath))
            queries = self._parse_queries(q_var.get())
            if queries:
                query_rows.append((name, queries, gpath))

        if not hit_dbs:
            self._set_output("No hit databases selected.")
            return None

        selected_dbtypes: dict[str, list[str]] = {}
        for name, _hdr, _sfx, _n, gpath in hit_dbs:
            selected_dbtypes.setdefault(dbtype_for(gpath), []).append(name)
        if len(selected_dbtypes) > 1:
            labels = {"nucl": "nucleotide", "prot": "protein"}
            lines = ["Error: Cannot mix nucleotide and protein databases as hit selections."]
            for dbtype, names in sorted(selected_dbtypes.items()):
                lines.append(f"  {labels.get(dbtype, dbtype)}: {', '.join(names)}")
            self._set_output("\n".join(lines))
            return None

        # Validate query IDs exist in their genome files
        errors = []
        for name, queries, gpath in query_rows:
            ids = load_fasta_ids(gpath)
            for q in queries:
                if q not in ids:
                    errors.append(f"  {q}  not found in {name}")
        if errors:
            self._set_output("Query ID errors:\n" + "\n".join(errors))
            return None

        # Check at least one query is provided
        if not query_rows:
            self._set_output("Error: At least one query (-q) must be specified.")
            return None

        # Collect -q / -qdbs pairs
        all_queries = []
        all_qdbs = []
        for name, queries, _gpath in query_rows:
            for q in queries:
                all_queries.append(q)
                all_qdbs.append(name)

        dbs = [s[0] for s in hit_dbs]
        hdrs = [s[1] for s in hit_dbs]
        sfxs = [s[2] or "none" for s in hit_dbs]
        ns = [s[3] for s in hit_dbs]

        has_suffix = any(s != "none" for s in sfxs)

        parts = []
        parts.extend(["-q"] + all_queries)
        parts.extend(["-qdbs"] + all_qdbs)
        parts.extend(["-dbs"] + dbs)
        parts.extend(["-hdr"] + hdrs)
        if has_suffix:
            parts.extend(["-hdr_sfx"] + sfxs)
        parts.extend(["-n"] + ns)

        aligner = self.aligner_var.get()
        if aligner != "clustalo":
            parts.extend(["--aligner", aligner])
        if aligner == "mafft" and self.mafft_mode_var.get() != "auto":
            parts.extend(["--mafft_mode", self.mafft_mode_var.get()])
        if self.tree_var.get() != "FastTree":
            parts.extend(["--tree_builder", self.tree_var.get()])
        if self.blast_type_var.get() != "tblastn":
            parts.extend(["--blast_type", self.blast_type_var.get()])

        threads = self.threads_var.get().strip()
        if threads and threads != str(max(1, os.cpu_count() // 2)):
            parts.extend(["--threads", threads])

        # Advanced options
        add_seqs = self.add_seqs_var.get().split()
        if add_seqs:
            parts.extend(["-add"] + add_seqs)

        add_dbs = self.add_dbs_var.get().split()
        if add_dbs:
            parts.extend(["-add_db"] + add_dbs)

        reroot = self.reroot_var.get().split()
        if len(reroot) > 1:
            self._set_output("Error: -a takes a single tip to root on, "
                             f"got {len(reroot)}: {' '.join(reroot)}")
            return None
        if reroot:
            parts.extend(["-a", reroot[0]])

        aa_start = self.aa_start_var.get().strip()
        aa_end = self.aa_end_var.get().strip()
        if aa_start and aa_end:
            parts.extend(["-aa", aa_start, aa_end])

        motifs = self.motif_var.get().split()
        if motifs:
            parts.extend(["--motif"] + motifs)
            if self.motif_syntax_var.get() != "regex":
                parts.extend(["--motif_syntax", self.motif_syntax_var.get()])
            if self.motif_overlap_var.get():
                parts.append("--motif_overlap")

        hmms = self.hmm_var.get().split()
        if hmms:
            parts.extend(["--hmm"] + hmms)

        return parts

    def _generate(self):
        self._save_settings()
        parts = self._build_args()
        if parts is None:
            return

        # Build display string
        display_parts = []
        i = 0
        while i < len(parts):
            if parts[i].startswith("-"):
                # Collect flag and its values
                flag = parts[i]
                vals = []
                i += 1
                while i < len(parts) and not parts[i].startswith("-"):
                    vals.append(parts[i])
                    i += 1
                if vals:
                    display_parts.append(flag + " " + " ".join(vals))
                else:
                    display_parts.append(flag)
            else:
                display_parts.append(parts[i])
                i += 1

        result = " ".join(display_parts)
        if self.prepend_var.get():
            result = "blast-align-tree " + result
        self._set_output(result)
        self.notebook.select(0)  # Switch to Command tab

    def _copy_to_clipboard(self):
        text = self.output.get("1.0", "end-1c")
        if text:
            self.root.clipboard_clear()
            self.root.clipboard_append(text)
            self.root.update()

    def _set_output(self, text: str):
        self.output.delete("1.0", "end")
        self.output.insert("1.0", text)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Tkinter GUI for building blast-align-tree commands."
    )
    parser.add_argument(
        "--default-n",
        type=positive_int,
        default=None,
        help=("default number of BLAST hits per selected genome "
              f"(default: {DEFAULT_N}, or the saved value from the last session)"),
    )
    parser.add_argument(
        "--select-genome",
        action="append",
        default=[],
        metavar="PATH",
        help=("preselect a genome by its path relative to genomes/; "
              "repeat this option to select more than one. Overrides the "
              "selection saved from the last session"),
    )
    parser.add_argument(
        "--no-restore",
        action="store_true",
        help=f"start from defaults, ignoring {SETTINGS_FILE.name}",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None):
    args = parse_args(argv)
    root = tk.Tk()
    GenomeSelectorApp(
        root,
        default_n=args.default_n,
        selected_genomes=args.select_genome,
        restore=not args.no_restore,
    )
    root.mainloop()


if __name__ == "__main__":
    main()

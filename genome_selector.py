#!/usr/bin/env python3
"""
BAT Genome Selector — tkinter GUI for building blast_align_tree.py commands.

Scans genomes/ for *.fa and *.fna files, auto-detects header tokens from the
first FASTA record, and generates copy-paste-ready -dbs / -hdr / -n arguments.
"""

import argparse
import os
import platform
import shutil
import random
import re
import subprocess
import tkinter as tk
from tkinter import ttk
from pathlib import Path


PROJ_DIR = Path(__file__).parent
GENOMES_DIR = PROJ_DIR / "genomes"
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
        with open(filepath) as f:
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
        with open(filepath) as f:
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
        with open(filepath) as f:
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
        with open(filepath) as f:
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
        with open(filepath) as f:
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
    def __init__(self, root: tk.Tk, default_n: int = DEFAULT_N):
        self.root = root
        root.title("BAT Genome Selector")
        root.minsize(1050, 600)

        self.default_n = default_n

        self._build_banner()
        self._build_genome_table()
        self._build_controls()
        self._build_options_panel()
        self._build_advanced_panel()
        self._build_action_buttons()
        self._build_output_panel()

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
        banner = (
            "  B.A.T.  BLAST - ALIGN - TREE\n"
            "  ---------------------------------------------------------------\n"
            "      .---------------- ATGCATGC  H. sapiens        /\\       /\\\n"
            "  .---+----.                                      __/  \\_____/  \\__\n"
            "  |        `------ ATGCGTGC  D. melanogaster    /  _   B.A.T.  _  \\\n"
            "  |                                             /__/ \\_________/ \\__\\\n"
            "  `---.---------- ATGCATCC  A. thaliana            \\___     ___/\n"
            "      `---------- ATGCGTCC  O. sativa                  \\___/\n"
            "\n"
            "  build the command, then run it in the environment you want"
        )
        tk.Label(self.root, text=banner, font=("Courier", 9), justify="left",
                 anchor="w", padx=12, pady=4).pack(fill="x")

        exts = "  ".join(sorted(FASTA_EXTENSIONS))
        summary = (
            f"Scanning genomes/ for FASTA files ({exts})\n"
            "Select genomes, configure headers, and generate a blast_align_tree.py command."
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

        # Sync header column widths with data rows
        self.inner_frame.bind("<Configure>", lambda e: (
            canvas.configure(scrollregion=canvas.bbox("all")),
            self._sync_columns(),
        ))

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
        ttk.Entry(row1, textvariable=self.add_dbs_var, width=30).pack(side="left", padx=(2, 0))
        ToolTip(row1.winfo_children()[-1], "Space-separated FASTA files for outgroup sequences")

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
        ttk.Checkbutton(gen_frame, text="Prepend  python blast_align_tree.py",
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
        self.notebook.pack(fill="both", expand=True, padx=8, pady=(0, 8))

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

        # Scrollable list
        list_frame = tk.Frame(runs_frame)
        list_frame.pack(fill="both", expand=True, padx=4, pady=(0, 4))

        self.runs_canvas = tk.Canvas(list_frame, highlightthickness=0)
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

    def _refresh_runs(self):
        """Populate the Recent Runs list."""
        for widget in self.runs_inner.winfo_children():
            widget.destroy()

        runs = scan_recent_runs()
        if not runs:
            tk.Label(self.runs_inner, text="No runs found.",
                     font=("TkDefaultFont", 9), fg="gray").grid(
                row=0, column=0, padx=8, pady=8)
            return

        # Column headers
        for col, text in enumerate(["Entry", "Timestamp", "Contents", ""]):
            tk.Label(self.runs_inner, text=text,
                     font=("TkDefaultFont", 9, "bold")).grid(
                row=0, column=col, padx=6, pady=(4, 2), sticky="w")

        for i, (entry, timestamp, path) in enumerate(runs, start=1):
            # Format timestamp: 20260131_0814 -> 2026-01-31 08:14
            display_ts = timestamp
            if len(timestamp) >= 13 and "_" in timestamp:
                try:
                    display_ts = (f"{timestamp[:4]}-{timestamp[4:6]}-{timestamp[6:8]}"
                                  f" {timestamp[9:11]}:{timestamp[11:13]}")
                except (IndexError, ValueError):
                    pass

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
        """Match header_frame column widths to inner_frame data columns."""
        for col in range(self.inner_frame.grid_size()[0]):
            if col == DIVIDER_COL:
                continue
            try:
                slaves = self.inner_frame.grid_slaves(row=0, column=col)
                if slaves:
                    w = slaves[0].winfo_width()
                    self.header_frame.columnconfigure(col, minsize=w)
            except (tk.TclError, IndexError):
                pass
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

        for i, gpath in enumerate(genome_files):
            tokens = detect_header_tokens(gpath)
            if "id" in tokens:
                tokens.remove("id")
            tokens.insert(0, "id")
            default_hdr = "id"

            chk_var = tk.BooleanVar(value=False)
            hdr_var = tk.StringVar(value=default_hdr)
            sfx_var = tk.StringVar(value="none")
            n_var = tk.StringVar(value=str(self.default_n))
            q_var = tk.StringVar(value="")
            example_header = get_example_header(gpath)
            preview_var = tk.StringVar(
                value=parse_header_token(example_header, default_hdr, sfx_var.get()))

            rel = gpath.relative_to(GENOMES_DIR).as_posix()
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
        _id_cache.clear()
        for widget in self.inner_frame.winfo_children():
            widget.destroy()
        self._populate_rows()
        self.root.update_idletasks()
        self._sync_columns()
        self._update_selected_display()

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
        """Build the argument list for blast_align_tree.py. Returns None on error."""
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
            result = "python blast_align_tree.py " + result
        self._set_output(result)
        self.notebook.select(0)  # Switch to Command tab
        print(result)

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
        description="Tkinter GUI for building blast_align_tree.py commands."
    )
    parser.add_argument(
        "--default-n",
        type=positive_int,
        default=DEFAULT_N,
        help=f"default number of BLAST hits per selected genome (default: {DEFAULT_N})",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None):
    args = parse_args(argv)
    root = tk.Tk()
    GenomeSelectorApp(root, default_n=args.default_n)
    root.mainloop()


if __name__ == "__main__":
    main()

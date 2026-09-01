# Unreleased — Explicit, auditable handling of identifier collisions and internal stop codons

> Draft. Version number not yet assigned; more features are still landing.
> Rename to `RELEASE_NOTES_v<x.y.z>.md` and bump `pyproject.toml` at release time.

Two places where the pipeline silently discarded or fabricated sequence are made
explicit, documented, and logged: which record wins an identifier collision, and
what happens to an in-frame stop codon.

## Part 1 — Identifier collisions

Header parsing (`-hdr gene:`, `-hdr locus=`, …) maps FASTA descriptions onto
short tree-tip identifiers, and that mapping is not one-to-one. Previously,
records that parsed to the same identifier were silently reduced to whichever
one happened to come first, discarding distinct isoforms, duplicated loci, or
records from different genomes sharing a gene symbol — with no record of what
was lost. This release makes every one of those decisions explicit, deterministic,
and logged.

The scale is not marginal: under `-hdr gene:`, 10,690 of the 48,321 records in
the bundled `TAIR10cds.fa` share a locus identifier with another record.

### Highlights

- **Deterministic isoform policy.** When several source records in one database
  claim the same identifier, the **longest amino-acid sequence** is retained;
  ties break on the lexicographically smallest source ID. This replaces
  first-in-sorted-glob order, which was effectively arbitrary. Note that the
  retained isoform is not necessarily the one passed as `-q`.
- **Genomes are never dropped.** When two *databases* use the same identifier,
  both records survive with a genome tag (`RPS5_TAIR10cds`, `RPS5_NbLab360`)
  instead of one silently replacing the other. Tags are applied only to
  contested identifiers, so ordinary tip labels are unchanged.
- **Per-query confirmation.** With a terminal attached, each query's collisions
  are shown for review after BLAST finishes but *before* results from multiple
  queries are merged — so genuine data loss is confirmed separately from the
  redundancy that merging queries is expected to produce. Answering `n` keeps
  every record, disambiguated as `AT2G13800__1`, `AT2G13800__2`, ….
- **Prompts time out after 10 seconds** and take their default, so a run left
  unattended finishes instead of blocking at the question. A default taken in
  silence is not passed off as agreement: the unanswered queries are named in
  the end-of-run summary and in a `# unanswered prompts (default taken):` line
  in `deduplication_log.tsv`.
- **New `--duplicates {ask,auto,fail}`.** Defaults to `ask`, falling back to
  `auto` when there is no terminal, so unattended and HPC runs never block.
  `fail` stops the run and lists the collisions.
- **`deduplication_log.tsv` beside the tree PDFs.** One row per affected record:
  stage, query, database, identifier, action (`kept` / `renamed` / `dropped`),
  source ID, amino-acid length, the original FASTA header, and the reason. A
  comment header records the run's `-hdr` / `-hdr_sfx` rule per database, so the
  file is self-describing.
- **End-of-run summary**, with collisions and cross-query overlap reported as
  two separate blocks so they are never conflated.
- **Accounting check.** The run aborts if the merged tree ends up with fewer
  tips than the resolution accounted for — a record cannot be lost outside the
  audited path.


### Compatibility

Runs with no identifier collisions are unaffected: with `-hdr id`, outputs are
byte-identical to v1.0.1 (verified end-to-end on `TAIR10cds.fa` — merged FASTA,
alignment, Newick tree, genome mapping, and per-genome hit files), and the only
new file is `deduplication_log.tsv`.

Runs that *did* have collisions may now produce a different representative
sequence per locus than before, because the previous winner was determined by
glob order rather than by any rule. Tip labels are unchanged unless an
identifier was contested across databases.

### Full changelog

- New `blast_align_tree/identifiers.py`: collision index, classification,
  resolution policy, and TSV ledger
- `cli.py`: identifier reconciliation stage between BLAST and the cross-query
  merge; parsed FASTAs and coding tables regenerated from the resolved map as a
  single enforcement point
- `cli.py`: `--duplicates` flag; end-of-run collision and query-overlap reports
- `cli.py`: `-add` sequences are checked against resolved identifiers and tagged
  rather than overwriting a BLAST hit
- `_parse_header_token` moved to `identifiers.py` as `parse_header_token`
- New `tests/test_identifiers.py` covering the policy

---

## Part 2 — Internal stop codons, reading frame, and strand

Previously, in-frame stop codons were **deleted** from every `tblastn` hit and
the flanking sequence joined. That produced a protein the input genome does not
encode, and — because deletion joins the reading frame on either side of the stop
— a pseudogene came out of the pipeline looking intact. The transformation was
also applied to the files labelled as nucleotide output, so no untransformed copy
of any hit survived the run.

### Highlights

- **New `--internal-stops {truncate,readthrough,excise}`, defaulting to
  `truncate`.** Translation now ends at the first in-frame stop, so the reported
  protein is the truncated product the locus encodes. `readthrough` keeps every
  codon in register and writes each stop as `X`, for domain-architecture work
  that needs downstream content. `excise` is the v1.0 behaviour, retained only to
  reproduce older runs.
- **Nucleotide outputs are now the untransformed records.** `all_hits.nt.fa` and
  `<entry>.nt.parse.merged.fa` are exactly what `blastdbcmd` returned, under
  every policy. Verified byte-for-byte against the BLAST database.
- **Reading frame and strand are diagnosed, never inferred.** Forward frame +1 of
  the record as stored is still always what gets translated. All six frames are
  translated for diagnosis only, and records that read substantially further in
  another frame are flagged `possible_utr` or `possible_wrong_strand`. Flags are
  warnings, never corrections — nothing is silently re-framed.
- **`translation_report.tsv` beside the tree PDFs.** One row per sequence: number
  of internal stops, position of the first, codons of reading frame remaining
  after it, longest stop-free run downstream, frame used, best-reading frame, and
  flags. A comment header records the policy in force and defines every flag.
- **Amino-acid coordinates are back in register with the genome.** Excision
  shifted every downstream residue by one position per excised stop, displacing
  HMM and motif coordinates. Under truncation, residue *p* maps to nucleotides
  `3p−2..3p` of the retrieved record.
- **The stop policy cannot change which isoform survives a collision.** Isoform
  ranking uses a stop-independent length, so a pseudogene is not silently
  superseded by an intact isoform just because truncation made it shorter.
- **Query translations are diagnosed but not altered**, since changing them would
  change the search itself.

### Compatibility

For CDS databases the new default is a no-op. Re-running the tutorial query
against `TAIR10cds.fa`, `Pvul218cds.fa` and `Vung469cds.fa` (45 hits) reports
zero internal stops, and the merged FASTA, alignment, Newick tree, and all three
per-genome amino-acid files are byte-identical between `truncate` and `excise`.

Output changes only for records that contain an internal stop, where the protein
is now truncated rather than fused across the stop. Such records are always
flagged in `translation_report.tsv`: verified on a synthetic transcript database
that no sequence changed without carrying a flag.

**Under truncation, domains downstream of an internal stop are no longer
detected**, so a degraded locus renders as a short domain architecture rather
than a continuous one. Use `--internal-stops readthrough` where downstream
content is needed.

### Full changelog

- New `blast_align_tree/translation.py`: policies, six-frame diagnosis, flag
  vocabulary, per-record stats, TSV report, end-of-run summary
- `cli.py`: `--internal-stops` flag; `blastdbcmd` output no longer transformed;
  translation stats collected per BLAST job as sidecars and merged into the run
  report; `-add` sequences follow the run policy; query translations diagnosed
- `cli.py`: `_remove_stop_codons()` removed; the legacy excision now lives in
  `translation.excise_stop_codons()` as the single definition of that behaviour
- `identifiers.build_index()`: optional `ranking_lengths` so isoform ranking is
  independent of the stop policy
- New `tests/test_translation.py` covering the policy, the flags, the register
  property, and legacy equivalence


## Part 3 — Runs record the command that made them

Until now an archived run in `ENTRY/runs/<timestamp>/` held its PDFs, trees and
logs but no record of how it was produced. `-n`, the aligner, the tree builder
and `-add` leave no trace in any output file, so a run from three weeks ago
could not be reproduced or even described without the shell history that
launched it.

- **New `run_command.txt` in every run.** Written from `sys.argv` at archive
  time — what was actually typed, not a re-rendering of argparse defaults — and
  reported in the end-of-run summary beside the alignment, tree and PDFs.
- **The BAT Genome Selector's Recent Runs tab reads it.** Selecting a run shows
  its queries, its databases with the `-n` used for each, its BLAST type,
  aligner and tree builder (marking the ones left at their defaults), plus
  extras such as `-add`, `--motif`, `--hmm` and `--datasets`.
- **Two copy buttons per run row.** `Re-run` copies the original command;
  `Re-draw` copies the `Rscript … visualize_tree.r … -n <NODE>` hint the
  pipeline prints when it finishes, carrying `--datasets` over from the
  original run. `Re-run` is disabled, with a tooltip, for runs that have no
  recorded command.
- **The tab also summarizes the run itself**: tips in the combined tree broken
  down per genome, de-duplication and translation log counts, output PDFs, and
  the basenames of any re-draws already made.

Runs archived before this release have no `run_command.txt`; the tab says so
and disables that copy rather than reconstructing a partial command from the
logs, which would silently omit `-n` and `-add` and so produce a different tree.

### Full changelog

- `cli.py`: `write_run_command()`; `RUN_COMMAND_NAME`; the file is written
  beside `deduplication_log.tsv` and listed in the end-of-run summary
- `genome_selector.py`: `read_run_command()`, `command_flags()`,
  `describe_command()`, `tree_stats()`, `log_stats()`, `redraw_names()`,
  `redraw_command()`, `format_timestamp()`; Recent Runs gains a details panel
  with per-row "Details" selection
- New `tests/test_run_details.py` covering the round trip and each summary

## Part 4 — Rooting the tree on an outgroup from the pipeline

`visualize_tree.r` has always had `-a/--reroot`, but the pipeline called it with
only `--entry`/`--write`/`--datasets`, so a run that pulled an outgroup in with
`-add`/`-add_db` still drew its first-pass PDFs unrooted; rooting them meant a
second, manual `Rscript` pass. `blast-align-tree` now takes `-a/--reroot ID` and
passes it straight through, so the PDFs a run produces come out already rooted.

### Highlights

- **`-a/--reroot ID` on the pipeline.** One tip, forwarded to
  `visualize_tree.r --reroot`. Pairs with `-add`/`-add_db`, which is what puts
  the outgroup in the tree in the first place.
- **A mistyped ID no longer wastes the run.** The R script matches tip labels
  exactly and errors out on no match — at the last step of a pipeline that may
  have run for an hour. The ID is now checked against the tree first; if it is
  not a tip, the run warns, draws the tree unrooted, and finishes, so the
  alignment, tree, and PDFs all survive.
- **The warning names the likely fix.** Tip labels are whatever `-hdr` parsing
  leaves behind, so the usual mistake is an isoform suffix the labels do not
  carry: `-a AT2G38240.1` answers with `Did you mean: AT2G38240?`.
- **Carried into the redraw hint and the run summary.** The `Rscript …` command
  printed at the end of the run repeats `-a`, and the Recent Runs details panel
  lists the reroot alongside the other extras.
- **In the selector too.** `bat-genome-selector` gains an `-a (reroot on)`
  field on the Advanced panel's outgroup row, saved with the rest of the
  settings. More than one tip is refused rather than silently truncated.

### Full changelog

- `cli.py`: `-a/--reroot` argument; `read_tip_labels()`, `resolve_reroot()`;
  `visualize_tree()` gains a `reroot` parameter; the end-of-run redraw hint
  repeats `-a`
- `genome_selector.py`: `reroot_var` entry on the Advanced panel, persisted
  via `GLOBAL_SETTING_VARS`, emitted by `_build_args()`; `describe_command()`
  reports `reroot` in Extras, and now matches long flag forms (`--add_seqs`,
  `--motifs`, `--hmms`, `--slice`) as well as short ones
- New `tests/test_reroot.py` covering the tip lookup, the near-match warning,
  and the flag's path into the Rscript command

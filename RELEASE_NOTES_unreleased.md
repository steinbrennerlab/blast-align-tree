# Unreleased — Explicit, auditable handling of identifier collisions

> Draft. Version number not yet assigned; more features are still landing.
> Rename to `RELEASE_NOTES_v<x.y.z>.md` and bump `pyproject.toml` at release time.

Header parsing (`-hdr gene:`, `-hdr locus=`, …) maps FASTA descriptions onto
short tree-tip identifiers, and that mapping is not one-to-one. Previously,
records that parsed to the same identifier were silently reduced to whichever
one happened to come first, discarding distinct isoforms, duplicated loci, or
records from different genomes sharing a gene symbol — with no record of what
was lost. This release makes every one of those decisions explicit, deterministic,
and logged.

The scale is not marginal: under `-hdr gene:`, 10,690 of the 48,321 records in
the bundled `TAIR10cds.fa` share a locus identifier with another record.

## Highlights

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

## Compatibility

Runs with no identifier collisions are unaffected: with `-hdr id`, outputs are
byte-identical to v1.0.1 (verified end-to-end on `TAIR10cds.fa` — merged FASTA,
alignment, Newick tree, genome mapping, and per-genome hit files), and the only
new file is `deduplication_log.tsv`.

Runs that *did* have collisions may now produce a different representative
sequence per locus than before, because the previous winner was determined by
glob order rather than by any rule. Tip labels are unchanged unless an
identifier was contested across databases.

## Full changelog

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

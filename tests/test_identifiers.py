"""Policy tests for identifier collision handling.

These exercise the rules the reviewer asked to be made explicit: which record
survives a collision, that the choice is deterministic, that nothing disappears
without a log row, and that benign overlap between queries is not mistaken for
data loss. No BLAST run is required.
"""

import pytest

from blast_align_tree.identifiers import (
    CROSS_DATABASE,
    CROSS_QUERY,
    ISOFORM,
    RecordRef,
    build_log_rows,
    collisions_by_query,
    detect,
    genome_tag,
    parse_header_token,
    query_overlaps,
    resolve,
    summarize,
)

AT = "TAIR10cds.fa"
NB = "NbLab360.v103.gff3.CDS.fasta"


def ref(query, db, source_id, aa_len, token):
    return RecordRef(query=query, db=db, source_id=source_id,
                     description=f"{source_id} gene:{token}", aa_len=aa_len, token=token)


def resolve_index(index, declined=()):
    collisions = detect(index)
    return collisions, resolve(index, collisions, declined)


# -----------------------
# Header token parsing
# -----------------------

def test_parse_header_token_matches_real_tair_header():
    desc = ("AT2G27490.4 cds chromosome:TAIR10:2:11748087:11749153:-1 gene:AT2G27490 "
            "gene_biotype:protein_coding gene_symbol:COAE")
    assert parse_header_token(desc, "gene:", "AT2G27490.4") == "AT2G27490"


def test_parse_header_token_falls_back_to_id_when_rule_absent():
    assert parse_header_token("Niben101Ctg00054g00001.1", "gene:", "Niben101Ctg00054g00001.1") \
        == "Niben101Ctg00054g00001.1"


def test_parse_header_token_strips_suffix():
    assert parse_header_token("x gene:AT1G01010.1 y", "gene:", "x", suffix=".1") == "AT1G01010"


# -----------------------
# Isoform collisions
# -----------------------

def test_longest_isoform_wins():
    index = [
        ref("RLK1", AT, "AT2G27490.1", 412, "AT2G27490"),
        ref("RLK1", AT, "AT2G27490.4", 388, "AT2G27490"),
    ]
    collisions, res = resolve_index(index)

    assert [c.kind for c in collisions] == [ISOFORM]
    assert res.final_ids == {(AT, "AT2G27490.1"): "AT2G27490"}
    assert res.dropped == {(AT, "AT2G27490.4")}


@pytest.mark.parametrize("order", [(0, 1), (1, 0)])
def test_winner_is_independent_of_input_order(order):
    records = [
        ref("RLK1", AT, "AT2G27490.4", 388, "AT2G27490"),
        ref("RLK1", AT, "AT2G27490.1", 412, "AT2G27490"),
    ]
    _, res = resolve_index([records[i] for i in order])
    assert res.final_ids == {(AT, "AT2G27490.1"): "AT2G27490"}


def test_equal_lengths_break_on_smallest_source_id():
    index = [
        ref("RLK1", AT, "AT2G27490.4", 400, "AT2G27490"),
        ref("RLK1", AT, "AT2G27490.1", 400, "AT2G27490"),
    ]
    _, res = resolve_index(index)
    assert res.final_ids == {(AT, "AT2G27490.1"): "AT2G27490"}


def test_three_isoforms_keep_exactly_one():
    index = [ref("RLK1", AT, f"AT2G27490.{i}", 300 + i, "AT2G27490") for i in (1, 2, 3)]
    _, res = resolve_index(index)
    assert res.final_ids == {(AT, "AT2G27490.3"): "AT2G27490"}
    assert len(res.dropped) == 2


# -----------------------
# Cross-query collisions
# -----------------------

def test_different_queries_hitting_different_isoforms_is_cross_query():
    index = [
        ref("RLK1", AT, "AT2G27490.1", 412, "AT2G27490"),
        ref("RLK2", AT, "AT2G27490.4", 388, "AT2G27490"),
    ]
    collisions, res = resolve_index(index)

    assert [c.kind for c in collisions] == [CROSS_QUERY]
    # Same rule as the within-query case, but no query is asked to confirm it.
    assert collisions_by_query(collisions) == {}
    assert res.dropped == {(AT, "AT2G27490.4")}


def test_same_record_hit_by_two_queries_is_overlap_not_collision():
    index = [
        ref("RLK1", AT, "AT2G27490.1", 412, "AT2G27490"),
        ref("RLK2", AT, "AT2G27490.1", 412, "AT2G27490"),
    ]
    collisions, res = resolve_index(index)

    assert collisions == []
    assert res.dropped == set()

    overlaps = query_overlaps(index)
    assert len(overlaps) == 1
    assert (overlaps[0].shared, overlaps[0].total_a, overlaps[0].total_b) == (1, 1, 1)


def test_query_needing_confirmation_is_the_one_that_self_collides():
    index = [
        ref("RLK1", AT, "AT2G27490.1", 412, "AT2G27490"),
        ref("RLK1", AT, "AT2G27490.4", 388, "AT2G27490"),
        ref("RLK2", AT, "AT2G27490.1", 412, "AT2G27490"),
    ]
    collisions, _ = resolve_index(index)
    assert sorted(collisions_by_query(collisions)) == ["RLK1"]


# -----------------------
# Cross-database collisions
# -----------------------

def test_shared_symbol_across_genomes_keeps_both_with_tags():
    index = [
        ref("RLK1", AT, "AT1G01010.1", 412, "RPS5"),
        ref("RLK1", NB, "NbL00g00010.1", 388, "RPS5"),
    ]
    collisions, res = resolve_index(index)

    assert [c.kind for c in collisions] == [CROSS_DATABASE]
    assert res.dropped == set()
    assert res.final_ids == {
        (AT, "AT1G01010.1"): "RPS5_TAIR10cds",
        (NB, "NbL00g00010.1"): "RPS5_NbLab360",
    }


def test_cross_database_tag_applies_only_to_colliding_identifiers():
    index = [
        ref("RLK1", AT, "AT1G01010.1", 412, "RPS5"),
        ref("RLK1", NB, "NbL00g00010.1", 388, "RPS5"),
        ref("RLK1", AT, "AT3G05780.1", 500, "LON3"),
    ]
    _, res = resolve_index(index)
    assert res.final_ids[(AT, "AT3G05780.1")] == "LON3"


def test_isoform_resolution_runs_before_genome_tagging():
    index = [
        ref("RLK1", AT, "AT1G01010.1", 412, "RPS5"),
        ref("RLK1", AT, "AT1G01010.2", 300, "RPS5"),
        ref("RLK1", NB, "NbL00g00010.1", 388, "RPS5"),
    ]
    _, res = resolve_index(index)

    assert res.dropped == {(AT, "AT1G01010.2")}
    assert res.final_ids == {
        (AT, "AT1G01010.1"): "RPS5_TAIR10cds",
        (NB, "NbL00g00010.1"): "RPS5_NbLab360",
    }


def test_genome_tag_shortens_versioned_filenames():
    assert genome_tag(NB) == "NbLab360"
    assert genome_tag("TAIR10cds.fa") == "TAIR10cds"
    assert genome_tag("animals/Pvul218cds.fa") == "Pvul218cds"


# -----------------------
# Declining the prompt
# -----------------------

def test_decline_keeps_every_record_with_rank_suffixes():
    index = [
        ref("RLK1", AT, "AT2G27490.1", 412, "AT2G27490"),
        ref("RLK1", AT, "AT2G27490.4", 388, "AT2G27490"),
    ]
    _, res = resolve_index(index, declined=[(AT, "AT2G27490")])

    assert res.dropped == set()
    assert res.final_ids == {
        (AT, "AT2G27490.1"): "AT2G27490__1",   # longest ranks first
        (AT, "AT2G27490.4"): "AT2G27490__2",
    }


def test_decline_combines_with_genome_tag():
    index = [
        ref("RLK1", AT, "AT1G01010.1", 412, "RPS5"),
        ref("RLK1", AT, "AT1G01010.2", 300, "RPS5"),
        ref("RLK1", NB, "NbL00g00010.1", 388, "RPS5"),
    ]
    _, res = resolve_index(index, declined=[(AT, "RPS5")])

    assert res.dropped == set()
    assert res.final_ids == {
        (AT, "AT1G01010.1"): "RPS5_TAIR10cds__1",
        (AT, "AT1G01010.2"): "RPS5_TAIR10cds__2",
        (NB, "NbL00g00010.1"): "RPS5_NbLab360",
    }


# -----------------------
# No record disappears silently
# -----------------------

def _distinct_records(index):
    return {(r.db, r.source_id) for r in index}


@pytest.mark.parametrize("declined", [(), [(AT, "AT2G27490")]])
def test_every_record_is_either_kept_or_logged_as_dropped(declined):
    index = [
        ref("RLK1", AT, "AT2G27490.1", 412, "AT2G27490"),
        ref("RLK1", AT, "AT2G27490.4", 388, "AT2G27490"),
        ref("RLK2", AT, "AT2G27490.1", 412, "AT2G27490"),
        ref("RLK1", NB, "NbL00g00010.1", 388, "AT2G27490"),
        ref("RLK2", NB, "NbL00g00020.1", 250, "OTHER"),
    ]
    collisions, res = resolve_index(index, declined=declined)

    kept = set(res.final_ids)
    assert kept | res.dropped == _distinct_records(index)
    assert kept & res.dropped == set()

    # Identifiers handed to the aligner must be unique, or the tree silently
    # merges tips again.
    assert len(set(res.final_ids.values())) == len(res.final_ids)

    rows = build_log_rows(index, res)
    for db, source_id in res.dropped:
        matching = [r for r in rows if r["source_id"] == source_id
                    and r["database"] == db and r["action"] == "dropped"]
        assert matching, f"{source_id} vanished without a log row"
        assert "superseded by" in matching[0]["reason"]


def test_overlap_rows_name_the_query_that_contributed_the_record():
    index = [
        ref("RLK1", AT, "AT3G05780.1", 500, "LON3"),
        ref("RLK2", AT, "AT3G05780.1", 500, "LON3"),
    ]
    _, res = resolve_index(index)
    rows = build_log_rows(index, res)

    overlap = [r for r in rows if r["stage"] == "merge_overlap"]
    assert len(overlap) == 1
    assert overlap[0]["query"] == "RLK2"
    assert "'RLK1'" in overlap[0]["reason"]
    assert "no data lost" in overlap[0]["reason"]


def test_summary_counts_each_collision_kind():
    index = [
        ref("RLK1", AT, "AT2G27490.1", 412, "AT2G27490"),
        ref("RLK1", AT, "AT2G27490.4", 388, "AT2G27490"),
        ref("RLK1", AT, "AT1G01010.1", 412, "RPS5"),
        ref("RLK1", NB, "NbL00g00010.1", 388, "RPS5"),
        ref("RLK1", NB, "NbL00g00030.1", 200, "SHARED"),
        ref("RLK2", NB, "NbL00g00040.1", 100, "SHARED"),
    ]
    collisions, res = resolve_index(index)
    counts = summarize(collisions, res)

    assert counts[ISOFORM] == 1
    assert counts[CROSS_QUERY] == 1
    assert counts[CROSS_DATABASE] == 1
    assert counts["declined"] == 0
    assert counts["records_dropped"] == 2


# -----------------------
# Inert when there are no collisions
# -----------------------

def test_unique_identifiers_are_left_completely_alone():
    index = [
        ref("RLK1", AT, "AT3G05780.1", 500, "AT3G05780"),
        ref("RLK1", AT, "AT1G24405.1", 239, "AT1G24405"),
        ref("RLK1", NB, "NbL00g00010.1", 388, "NbL00g00010.1"),
    ]
    collisions, res = resolve_index(index)

    assert collisions == []
    assert res.dropped == set()
    assert all(final == r.token for r, final in
               ((r, res.final_ids[(r.db, r.source_id)]) for r in index))
    assert build_log_rows(index, res) == []

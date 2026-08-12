"""Policy tests for internal stop codons and reading frame.

These exercise the rules the reviewer asked to be made explicit: that a stop
codon truncates the reported protein rather than being deleted from it, that the
pipeline says how much reading frame remained after it, that the frame actually
used is never inferred from the data, and that a record whose frame looks wrong
is flagged rather than quietly re-framed. No BLAST run is required.
"""

import random

import pytest

from Bio.Seq import Seq

from blast_align_tree import identifiers
from blast_align_tree.translation import (
    EMPTY_TRANSLATION,
    EXCISE,
    INTERNAL_STOP,
    LEN_NOT_MULTIPLE_OF_3,
    NO_START_CODON,
    NO_TERMINAL_STOP,
    ORF_MOSTLY_DOWNSTREAM,
    POSSIBLE_UTR,
    POSSIBLE_WRONG_STRAND,
    POLICIES,
    READTHROUGH,
    TRUNCATE,
    USED_FRAME,
    analyze,
    excise_stop_codons,
    translate_sequence,
)


# Sense codons only, so a constructed CDS has no accidental internal stop.
SENSE_CODONS = ["GCT", "CGT", "AAC", "GAT", "TGC", "CAA", "GAA", "GGT", "CAT",
                "ATC", "TTA", "AAA", "ATG", "TTC", "CCA", "TCT", "ACT", "TGG",
                "TAC", "GTT"]


def make_cds(n_codons=300, seed=7):
    """A well-formed CDS: ATG, sense codons, terminal stop, nothing internal.

    Realistic length matters. A wrong reading frame only becomes recognisable
    once it is long enough to accumulate stop codons, so toy 5-codon fixtures
    would exercise a regime the check is not making claims about.
    """
    rng = random.Random(seed)
    return "ATG" + "".join(rng.choice(SENSE_CODONS) for _ in range(n_codons)) + "TGA"


def with_internal_stop(cds, codon_index, stop="TAA"):
    codons = [cds[i:i + 3] for i in range(0, len(cds), 3)]
    codons[codon_index] = stop
    return "".join(codons)


def legacy_remove_stop_codons(seq):
    """
    The v1.0 excision, copied verbatim from the pre-change cli.py so the legacy
    policy is checked against an independent copy rather than against itself.
    """
    stops = {"TAG", "TGA", "TAA"}
    seq_list = list(str(seq))
    i = 0
    while i + 2 < len(seq_list):
        codon = "".join(seq_list[i:i + 3]).upper()
        if codon in stops:
            del seq_list[i:i + 3]
        else:
            i += 3
    return "".join(seq_list)


def legacy_translate(seq):
    """v1.0 translation: trim to whole codons, then translate."""
    seq = seq[:len(seq) - (len(seq) % 3)]
    return str(Seq(seq).translate())


# -----------------------
# A well-formed CDS is unaffected
# -----------------------

def test_clean_cds_translation_is_unchanged_from_v1():
    cds = make_cds()
    expected = legacy_translate(legacy_remove_stop_codons(cds))
    assert translate_sequence(cds, TRUNCATE) == expected


def test_clean_cds_raises_no_flags():
    assert analyze(make_cds()).flags == ()


def test_clean_cds_drops_only_the_terminal_stop():
    cds = make_cds(n_codons=10)
    aa = translate_sequence(cds, TRUNCATE)
    assert aa.startswith("M")
    assert "*" not in aa
    assert len(aa) == 11  # ATG + 10 codons; terminal stop not a residue


# -----------------------
# Internal stops truncate, and the remaining frame is reported
# -----------------------

def test_internal_stop_truncates_rather_than_splicing():
    cds = make_cds()
    pseudo = with_internal_stop(cds, 200)
    aa = translate_sequence(pseudo, TRUNCATE)

    assert len(aa) == 200
    assert "*" not in aa
    # The residues that survive are exactly the ones the genome encodes there.
    assert aa == translate_sequence(cds, TRUNCATE)[:200]
    # And nothing downstream was pulled forward, which is what excision did.
    assert not translate_sequence(cds, TRUNCATE).startswith(aa + "?")


def test_reports_how_much_reading_frame_remained():
    cds = make_cds()                     # 302 codons: ATG + 300 + stop
    pseudo = with_internal_stop(cds, 200)
    stat = analyze(pseudo)

    assert stat.n_internal_stops == 1
    assert stat.first_stop_aa_pos == 201          # 1-based
    assert INTERNAL_STOP in stat.flags
    # 302 codons total, 201 consumed through the stop itself.
    assert stat.aa_after_first_stop == 101
    # The stop-free run after it: up to the terminal stop.
    assert stat.longest_downstream_orf_aa == 100


def test_second_stop_shortens_the_longest_downstream_orf():
    cds = make_cds()
    single = analyze(with_internal_stop(cds, 200))
    pseudo = with_internal_stop(with_internal_stop(cds, 200), 290)
    stat = analyze(pseudo)

    assert stat.n_internal_stops == 2
    assert stat.first_stop_aa_pos == 201
    assert translate_sequence(pseudo, TRUNCATE) == translate_sequence(cds, TRUNCATE)[:200]
    # Total frame remaining is unchanged, but it is now broken into two runs:
    # codons 202-290 and 292-301. The longest is the metric that matters for
    # judging whether a downstream ORF survives.
    assert stat.aa_after_first_stop == single.aa_after_first_stop == 101
    assert single.longest_downstream_orf_aa == 100
    assert stat.longest_downstream_orf_aa == 89


def test_no_terminal_stop_is_flagged():
    cds = make_cds()[:-3]
    assert NO_TERMINAL_STOP in analyze(cds).flags


# -----------------------
# Frame and strand are diagnosed, never corrected
# -----------------------

def test_out_of_frame_5_prime_utr_is_flagged_as_possible_utr():
    rng = random.Random(11)
    utr = "".join(rng.choice("ACGT") for _ in range(64))   # 64 % 3 != 0
    stat = analyze(utr + make_cds())

    assert POSSIBLE_UTR in stat.flags
    assert stat.best_frame == "+2"
    assert any(f.startswith("frame_mismatch") for f in stat.flags)


def test_flagged_record_is_still_translated_in_frame_plus_one():
    """The flag must not become a silent correction."""
    rng = random.Random(11)
    utr = "".join(rng.choice("ACGT") for _ in range(64))
    seq = utr + make_cds()
    stat = analyze(seq)

    assert stat.best_frame == "+2"          # the frame it *would* read better in
    assert USED_FRAME == "+1"
    # What is emitted is frame +1 regardless.
    assert translate_sequence(seq, TRUNCATE) == str(Seq(seq[:len(seq) // 3 * 3]).translate()).split("*", 1)[0]


def test_in_frame_5_prime_utr_does_not_trip_the_frame_flag():
    """A UTR whose length is a multiple of 3 leaves the reading frame intact."""
    rng = random.Random(13)
    utr = "".join(rng.choice(SENSE_CODONS) for _ in range(21))   # 63 nt
    stat = analyze(utr + make_cds())

    assert POSSIBLE_UTR not in stat.flags
    assert stat.best_frame == "+1"
    assert NO_START_CODON in stat.flags       # still visibly not a bare CDS


def test_record_stored_on_the_opposite_strand_is_flagged():
    rev = str(Seq(make_cds()).reverse_complement())
    stat = analyze(rev)

    assert POSSIBLE_WRONG_STRAND in stat.flags
    assert stat.best_frame == "-1"


def test_in_frame_utr_carrying_a_stop_is_caught_by_the_orf_comparison():
    """
    The case the frame check cannot see: an in-frame 5' UTR that happens to
    contain a stop codon. Frame +1 still holds the longest ORF, so no frame
    mismatch is reported - but the protein truncates after a handful of residues
    while the real coding sequence sits downstream.
    """
    cds = make_cds()
    seq = "ATG" + "TAA" + "GCT" * 5 + cds       # stop at codon 2 of the record
    stat = analyze(seq)

    assert ORF_MOSTLY_DOWNSTREAM in stat.flags
    assert not stat.frame_suspect                     # frame +1 is genuinely best
    assert stat.first_stop_aa_pos == 2
    assert len(translate_sequence(seq, TRUNCATE)) == 1
    assert stat.longest_downstream_orf_aa > 300


def test_ordinary_pseudogene_is_not_flagged_orf_mostly_downstream():
    """A stop two-thirds through leaves less ahead than behind: no flag."""
    stat = analyze(with_internal_stop(make_cds(), 200))
    assert ORF_MOSTLY_DOWNSTREAM not in stat.flags


def test_pseudogene_is_not_mistaken_for_a_frame_problem():
    """One internal stop must not look like a wrong frame."""
    stat = analyze(with_internal_stop(make_cds(), 200))

    assert stat.flags == (INTERNAL_STOP,)
    assert stat.best_frame == "+1"
    assert not stat.frame_suspect


# -----------------------
# Degenerate input
# -----------------------

@pytest.mark.parametrize("seq", ["", "A", "AT"])
def test_short_records_survive(seq):
    stat = analyze(seq)
    assert EMPTY_TRANSLATION in stat.flags
    assert translate_sequence(seq, TRUNCATE) == ""


def test_length_not_a_multiple_of_three_is_flagged():
    stat = analyze(make_cds() + "AT")
    assert LEN_NOT_MULTIPLE_OF_3 in stat.flags


def test_unknown_policy_is_rejected():
    with pytest.raises(ValueError):
        translate_sequence(make_cds(), "delete-everything")


# -----------------------
# Readthrough and the legacy policy
# -----------------------

def test_readthrough_marks_stops_without_inventing_residues():
    cds = make_cds()
    pseudo = with_internal_stop(cds, 200)
    aa = translate_sequence(pseudo, READTHROUGH)

    assert aa[200] == "X"                       # the stop, marked as unknown
    assert len(aa) == 301                       # every codon kept, terminal stop dropped
    assert aa[:200] == translate_sequence(cds, TRUNCATE)[:200]
    assert aa[201:] == translate_sequence(cds, TRUNCATE)[201:]


def test_excise_reproduces_v1_output_exactly():
    pseudo = with_internal_stop(make_cds(), 200)
    assert excise_stop_codons(pseudo) == legacy_remove_stop_codons(pseudo)
    assert translate_sequence(pseudo, EXCISE) == \
        legacy_translate(legacy_remove_stop_codons(pseudo))


def test_excise_is_the_policy_that_fabricates_a_junction():
    """The behaviour the reviewer objected to, pinned so it cannot come back by default."""
    cds = make_cds()
    pseudo = with_internal_stop(cds, 200)
    excised = translate_sequence(pseudo, EXCISE)

    # Excision joins the two sides into one continuous protein with no gap.
    assert len(excised) == 300
    assert "*" not in excised
    # And it shifts everything downstream by exactly one residue per removed stop.
    assert excised[200:] == translate_sequence(cds, TRUNCATE)[201:]


# -----------------------
# Amino-acid to nucleotide register
# -----------------------

def _codon_at(nt, aa_pos):
    """The codon encoding 1-based residue `aa_pos` of the record as retrieved."""
    return str(Seq(nt[3 * (aa_pos - 1):3 * aa_pos]).translate())


def test_truncate_keeps_residues_in_register_with_the_nucleotide_record():
    pseudo = with_internal_stop(make_cds(), 200)
    aa = translate_sequence(pseudo, TRUNCATE)
    for pos in (1, 57, 150, len(aa)):
        assert aa[pos - 1] == _codon_at(pseudo, pos)


def test_excise_breaks_that_register_downstream_of_a_stop():
    """Why feature coordinates drifted under the legacy policy."""
    pseudo = with_internal_stop(make_cds(), 200)
    aa = translate_sequence(pseudo, EXCISE)

    assert aa[149] == _codon_at(pseudo, 150)        # upstream: still aligned
    assert aa[249] != _codon_at(pseudo, 250)        # downstream: off by the excised stop
    assert aa[249] == _codon_at(pseudo, 251)        # by exactly one residue


# -----------------------
# The stop policy must not change which isoform survives a collision
# -----------------------

def test_ranking_length_is_identical_under_every_policy():
    pseudo = with_internal_stop(make_cds(), 200)
    lengths = {p: analyze(pseudo, policy=p).aa_len_ranking for p in POLICIES}
    assert len(set(lengths.values())) == 1, lengths


def test_pseudogene_isoform_still_wins_its_collision_under_truncation():
    """
    A locus with two isoforms: the longer one carries an early stop. Ranking on
    the truncated length would silently discard it in favour of the shorter
    intact isoform - exactly the evidence a pseudogene claim rests on.
    """
    long_pseudo = with_internal_stop(make_cds(n_codons=300, seed=7), 20)
    short_intact = make_cds(n_codons=100, seed=9)

    def index_for(policy):
        return [
            identifiers.RecordRef(
                query="Q", db="db.fa", source_id=sid,
                description=f"{sid} gene:LOCUS",
                aa_len=analyze(seq, policy=policy).aa_len_ranking,
                token="LOCUS")
            for sid, seq in (("iso.1", long_pseudo), ("iso.2", short_intact))
        ]

    winners = set()
    for policy in POLICIES:
        index = index_for(policy)
        collisions = identifiers.detect(index)
        assert len(collisions) == 1
        winners.add(collisions[0].winner.source_id)

    assert winners == {"iso.1"}

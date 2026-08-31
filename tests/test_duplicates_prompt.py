"""Tests for the timed identifier-collision prompt.

A collision prompt used to block forever, so an unattended run stalled at the
question instead of finishing. It now takes its default answer after
CONFIRM_TIMEOUT_SECONDS -- and records that the default was taken rather than
agreed to, since 'ask' alone would overstate what happened.
"""

import threading
import time

import pytest

from blast_align_tree import cli
from blast_align_tree import identifiers
from blast_align_tree.identifiers import RecordRef, detect, resolve, write_log

AT = "TAIR10cds.fa"


class Stdin:
    """A stdin that hands over prepared lines, then blocks like a terminal."""

    def __init__(self, *lines, delay=0.1):
        self.lines = list(lines)
        self.delay = delay
        self.exhausted = threading.Event()

    def readline(self):
        # Answer only after the prompt has cleared pre-typed input, the way a
        # person replying to a question they can actually see does.
        time.sleep(self.delay)
        if self.lines:
            return self.lines.pop(0)
        self.exhausted.set()
        threading.Event().wait()   # a terminal with nobody typing at it

    def isatty(self):
        return True


class ClosedStdin(Stdin):
    def readline(self):
        time.sleep(self.delay)
        return ""                  # EOF, e.g. stdin closed mid-run


@pytest.fixture(autouse=True)
def fresh_reader(monkeypatch):
    """The stdin reader is process-wide; give each test its own."""
    monkeypatch.setattr(cli, "_stdin_reader", None)
    monkeypatch.setattr(cli, "_stdin_lines", cli.queue.Queue())
    yield


def test_silence_accepts_after_the_timeout(monkeypatch, capsys):
    monkeypatch.setattr(cli.sys, "stdin", Stdin())
    accepted, answered = cli._prompt_yes("Deduplicate?", timeout=0.2)
    assert (accepted, answered) == (True, False)
    assert "no answer in 0s: deduplicating" in capsys.readouterr().out


def test_an_answer_still_wins_before_the_timeout(monkeypatch):
    monkeypatch.setattr(cli.sys, "stdin", Stdin("n\n"))
    assert cli._prompt_yes("Deduplicate?", timeout=2) == (False, True)


def test_empty_line_accepts_as_the_default(monkeypatch):
    monkeypatch.setattr(cli.sys, "stdin", Stdin("\n"))
    assert cli._prompt_yes("Deduplicate?", timeout=2) == (True, True)


def test_eof_accepts_without_counting_as_an_answer(monkeypatch):
    monkeypatch.setattr(cli.sys, "stdin", ClosedStdin())
    assert cli._prompt_yes("Deduplicate?", timeout=2) == (True, False)


def test_prompt_shows_the_countdown_and_its_default(monkeypatch, capsys):
    monkeypatch.setattr(cli.sys, "stdin", Stdin("y\n"))
    cli._prompt_yes("Deduplicate?", timeout=10)
    assert "Deduplicate? [Y/n] (10s -> Y)" in capsys.readouterr().out


def test_a_stale_answer_is_not_applied_to_the_next_query(monkeypatch, capsys):
    """Typing after one prompt gave up must not silently answer the next one."""
    stdin = Stdin()
    monkeypatch.setattr(cli.sys, "stdin", stdin)
    assert cli._prompt_yes("Query 1?", timeout=0.2) == (True, False)

    stdin.exhausted.wait(timeout=2)
    cli._stdin_lines.put("n\n")               # arrives after the first prompt gave up
    assert cli._prompt_yes("Query 2?", timeout=0.2) == (True, False)


def test_unanswered_prompts_are_named_in_the_dedup_log(tmp_path):
    index = [
        RecordRef(query="RLK1", db=AT, source_id="AT2G27490.1",
                  description="AT2G27490.1 gene:AT2G27490", aa_len=412, token="AT2G27490"),
        RecordRef(query="RLK1", db=AT, source_id="AT2G27490.4",
                  description="AT2G27490.4 gene:AT2G27490", aa_len=388, token="AT2G27490"),
    ]
    collisions = detect(index)
    res = resolve(index, collisions, ())
    log = write_log(
        tmp_path / "deduplication_log.tsv",
        identifiers.build_log_rows(index, res),
        entry="RLK1", blast_type="tblastn",
        hdr_rules_by_db={AT: ("gene:", "")}, mode="ask",
        unanswered=["RLK1"],
    )
    header = log.read_text(encoding="utf-8").splitlines()
    assert "# duplicates mode: ask" in header
    assert "# unanswered prompts (default taken): RLK1" in header


def test_answered_runs_leave_no_such_line(tmp_path):
    log = write_log(
        tmp_path / "deduplication_log.tsv", [],
        entry="RLK1", blast_type="tblastn",
        hdr_rules_by_db={AT: ("gene:", "")}, mode="ask",
    )
    assert "unanswered" not in log.read_text(encoding="utf-8")


def test_timeout_follows_the_constant_not_a_bound_default(monkeypatch):
    """CONFIRM_TIMEOUT_SECONDS stays authoritative when it is changed."""
    monkeypatch.setattr(cli, "CONFIRM_TIMEOUT_SECONDS", 0.2)
    monkeypatch.setattr(cli.sys, "stdin", Stdin())
    assert cli._prompt_yes("Deduplicate?") == (True, False)

"""Tests for the align module."""

from align import (
    align,
    get_edits,
    local_align,
    edit_dist
)


def test_edits() -> None:
    """Testing the edits() function."""
    assert get_edits(
        "acca-aagt--a", "a-caaatgtcca"
    ) == ("accaaagta", "acaaatgtcca", "MDMMIMMMMIIM")
    assert get_edits(
        "acgttcga", "aaa---aa"
    ) == ("acgttcga", "aaaaa", "MMMDDDMM")


def test_align() -> None:
    """Testing the align() function."""
    assert align(
        "accaaagta",
        "acaaatgtcca",
        "MDMMIMMMMIIM"
    ) == ("acca-aagt--a", "a-caaatgtcca")
    assert align(
        "a", "", "D"
    ) == ("a", "-")


def test_local_align() -> None:
    """Testing the local_align() function."""
    assert local_align(
        "accaaagta",
        "gtacaaatgtcca",
        2,
        "MDMMIMMMMIIM"
    ) == ("acca-aagt--a", "a-caaatgtcca")
    assert local_align(
        "a", "", 0, "D"
    ) == ("a", "-")


def test_edit_dist() -> None:
    """Testing the edit_dist() function."""
    assert edit_dist("accaaagta", "cgacaaatgtcca", 2, "MDMMIMMMMIIM") == 5
    assert edit_dist("acgttcga", "aaaaa", 0, "MMMDDDMM") == 6

"""Tests for the cigar module."""

from cigar import (
    cigar_to_edits,
    edits_to_cigar
)


def test_cigar_to_edits() -> None:
    """Testing the cigar_to_edits() function."""
    assert cigar_to_edits("1M1D1I1M1I1D") == "MDIMID"
    assert cigar_to_edits("2M2D2I2M2I2D") == "MMDDIIMMIIDD"
    assert cigar_to_edits("1M2D3I2M1I2D") == "MDDIIIMMIDD"
    assert cigar_to_edits("") == ""


def test_edits_to_cigar() -> None:
    """Testing the edits_to_cigar() function."""
    assert edits_to_cigar("MDIMID") == "1M1D1I1M1I1D"
    assert edits_to_cigar("MMDDIIMMIIDD") == "2M2D2I2M2I2D"
    assert edits_to_cigar("MDDIIIMMIDD") == "1M2D3I2M1I2D"
    assert edits_to_cigar("") == ""

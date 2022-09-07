from typing import Tuple
"""A module for translating between alignments and edits sequences."""


def get_edits(p: str, q: str) -> Tuple[str, str, str]:
    """Extract the edit operations from a pairwise alignment.

    Args:
        p (str): The first row in the pairwise alignment.
        q (str): The second row in the pairwise alignment.

    Returns:
        str: The list of edit operations as a string.

    >>> get_edits('ACCACAGT-CATA', 'A-CAGAGTACAAA')
    ('ACCACAGTCATA', 'ACAGAGTACAAA', 'MDMMMMMMIMMMM')

    """
    assert len(p) == len(q)

    out_p = ''
    out_q = ''
    edits = ''

    for i in range(len(p)):
        cp = p[i]
        cq = q[i]
        if cp != '-' and cq != '-':
            out_p += cp
            out_q += cq
            edits += 'M'
        if cp == '-':
            out_q += cq
            edits += 'I'
        if cq == '-':
            out_p += cp
            edits += 'D'


    return out_p, out_q, edits


def local_align(p: str, x: str, i: int, edits: str) -> Tuple[str, str]:
    """Align two sequences from a sequence of edits.

    Args:
        p (str): The read string we have mapped against x
        x (str): The longer string we have mapped against
        i (int): The location where we have an approximative match
        edits (str): The list of edits to apply, given as a string

    Returns:
        tuple[str, str]: The two rows in the pairwise alignment

    >>> local_align("ACCACAGTCATA", "GTACAGAGTACAAA", 2, "MDMMMMMMIMMMM")
    ('ACCACAGT-CATA', 'A-CAGAGTACAAA')

    """
    return align(p, x[i:], edits)


def align(p: str, q: str, edits: str) -> Tuple[str, str]:
    """Align two sequences from a sequence of edits.

    Args:
        p (str): The first sequence to align.
        q (str): The second sequence to align
        edits (str): The list of edits to apply, given as a string

    Returns:
        tuple[str, str]: The two rows in the pairwise alignment

    >>> align("ACCACAGTCATA", "ACAGAGTACAAA", "MDMMMMMMIMMMM")
    ('ACCACAGT-CATA', 'A-CAGAGTACAAA')

    """
    out_p = ''
    out_q = ''
    i_p = 0
    i_q = 0

    for e in edits:
        if (e == 'M'):
            out_p += p[i_p]
            out_q += q[i_q]
            i_p += 1
            i_q += 1
        if (e == 'I'):
            out_p += '-'
            out_q += q[i_q]
            i_q += 1
        if (e == 'D'):
            out_p += p[i_p]
            out_q += '-'
            i_p += 1
    

    return out_p, out_q


def edit_dist(p: str, x: str, i: int, edits: str) -> int:
    """Get the distance between p and the string that starts at x[i:]
    using the edits.

    Args:
        p (str): The read string we have mapped against x
        x (str): The longer string we have mapped against
        i (int): The location where we have an approximative match
        edits (str): The list of edits to apply, given as a string

    Returns:
        int: The distance from p to x[i:?] described by edits

    >>> edit_dist("accaaagta", "cgacaaatgtcca", 2, "MDMMIMMMMIIM")
    5
    """
    aligned_p, aligned_q = local_align(p, x, i, edits)
    dist = 0
    
    for i in range (len(aligned_p)):
        c_p = aligned_p[i]
        c_q = aligned_q[i]
        if c_p == '-' or c_q == '-':
            dist += 1
        elif c_p != c_q:
            dist += 1
        

    return dist
# Copyright 2019
# Author: Fabio Gutmann <https://github.com/fabio-gut>


def idx_to_array_index(idx: int, idxpos0: int):
    """
    Converts a nucleotide index to an 0-included array index
    :param idx: The index to convert
    :param idxpos0: The start index
    :return: The index of the element in the array
    >>> idx_to_array_index(10, 5)
    5
    >>> idx_to_array_index(10, -200)
    209
    >>> idx_to_array_index(1, -1)
    1
    """
    return idx - idxpos0 - (1 if (idxpos0 < 0 < idx or idx < 0 < idxpos0) else 0)


def array_index_to_idx(i: int, idxpos0: int):
    """
    Converts a nucleotide index to an 0-included array index
    :param i: The array index
    :param idxpos0: The start index
    :return: The index of the element in the array
    >>> array_index_to_idx(5, 5)
    10
    >>> array_index_to_idx(209, -200)
    10
    >>> array_index_to_idx(1, -1)
    1
    """
    return idxpos0 + i + (1 if (idxpos0 < 0 < i or i < 0 < idxpos0) else 0)


def add_sub_idx(i1: int, change: int):
    """
    Adds or subtracts an index from another one and skips 0 if encountered.
    :param i1: The start index
    :param change: The change (+ or -)
    :return: The new index
    >>> add_sub_idx(-10, 10)
    1
    >>> add_sub_idx(10, -10)
    -1
    >>> add_sub_idx(-5, 10)
    6
    """
    return i1 + change + (1 if i1 < 0 < i1+change+1 else -1 if i1+change-1 < 0 < i1 else 0)

# Copyright 2019
# Author: Fabio Gutmann <https://github.com/fabio-gut>

from typing import Tuple, List
from copomus.indexing import idx_to_array_index


def _filter_bp_type(bps: List[Tuple[int, int]], q: str, t: str, qidxpos0: int, tidxpos0: int, bp_type='') -> \
        List[Tuple[int, int]]:
    """
    Filters out base pairs of a specified type
    :param bps: List of tuples in form (q_index, t_index) representing the base pairs
    :param q: The query sequence
    :param t: The target sequence
    :param qidxpos0: Starting index for the query
    :param tidxpos0: Starting index for the target
    :param bp_type: The type of base pairs to filter out, for example GU, CG, AU, GC, ...
    :return: List of tuples in form (q_index, t_index) representing the base pairs
    >>> q, t = 'GCUACGAUC', 'UUUGCGAGCAGCUAGG'
    >>> bps = [(1,1),(2,2),(6,13),(8,15),(9,16)]
    >>> _filter_bp_type(bps, q, t, 1, 1, 'GU')
    [(2, 2), (9, 16)]
    >>> _filter_bp_type(bps, q, t, 1, 1, 'GC')
    [(1, 1), (2, 2), (6, 13), (8, 15)]
    """
    new_bps = []

    for q_index, t_index in bps:
        q_array_index = idx_to_array_index(q_index, qidxpos0)
        t_array_index = idx_to_array_index(t_index, tidxpos0)
        if f'{q[q_array_index]}{t[t_array_index]}' in [bp_type, bp_type[::-1]]:
            continue
        new_bps.append((q_index, t_index))
    return new_bps


def _neighbors_in_mfe(q_index: int, t_index: int, bps: List[Tuple[int, int]]) -> int:
    """
    Counts how many neighboring base pairs are in the MFE region
    :param q_index: Base pair index on the query
    :param t_index: Base pair index on the target
    :param bps: List of tuples in form (q_index, t_index) representing the base pairs
    :return: Count of neighboring base pairs that are in the MFE region
    >>> _neighbors_in_mfe(3, 5, [(1, 7), (6, 9), (2, 6), (4, 4), (7, 4)])
    2
    >>> _neighbors_in_mfe(3, 5, [(1, 7), (6, 9), (4, 4), (7, 4)])
    1
    """
    count = 0
    for q_o, t_o in [(+1, -1), (-1, +1)]:  # get neighbors by offset
        if (q_index+q_o, t_index+t_o) in bps:
            count += 1
    return count


def _neighbors_can_pair(q_index: int, t_index: int, q: str, t: str, qidxpos0: int, tidxpos0: int) -> int:
    """
    Counts how many neighboring base pairs can pair
    :param q_index: Base pair index on the query
    :param t_index: Base pair index on the target
    :param q: The query sequence
    :param t: The target sequence
    :param qidxpos0: Starting index for the query
    :param tidxpos0: Starting index for the target
    :return: Count of neighboring base pairs that can pair
    >>> _neighbors_can_pair(1, 5, 'GCAUCGAUC', 'CGUACGAUCGAUCC', 1, 1)
    0
    >>> _neighbors_can_pair(3, 3, 'GCAUCGAUC', 'CGUACGAUCGAUCC', 1, 1)
    1
    >>> _neighbors_can_pair(6, 9, 'GCAUCGAUC', 'CGUACGAUCGAUCC', 1, 1)
    2
    """
    count = 0
    q_array_index = idx_to_array_index(q_index, qidxpos0)
    t_array_index = idx_to_array_index(t_index, tidxpos0)
    pairable_bps = ['GC', 'CG', 'AU', 'UA', 'GU', 'UG']
    for q_o, t_o in [(+1, -1), (-1, +1)]:  # get neighbors by offset
        if q_array_index+q_o in range(len(q)) and t_array_index+t_o in range(len(t)):
            if f'{q[q_array_index+q_o]}{t[t_array_index+t_o]}' in pairable_bps:
                count += 1
    return count


def filter_gu(bps: List[Tuple[int, int]], q: str, t: str, qidxpos0: int, tidxpos0: int) -> List[Tuple[int, int]]:
    """Filters any GU or UG base pair"""
    return _filter_bp_type(bps, q, t, qidxpos0, tidxpos0, bp_type='GU')


def filter_au(bps: List[Tuple[int, int]], q: str, t: str, qidxpos0: int, tidxpos0: int) -> List[Tuple[int, int]]:
    """Filters any AU or UA base pair"""
    return _filter_bp_type(bps, q, t, qidxpos0, tidxpos0, bp_type='AU')


def filter_gc(bps: List[Tuple[int, int]], q: str, t: str, qidxpos0: int, tidxpos0: int) -> List[Tuple[int, int]]:
    """Filters ay GC or CG base pair"""
    return _filter_bp_type(bps, q, t, qidxpos0, tidxpos0, bp_type='GC')


def filter_lp(bps: List[Tuple[int, int]], q: str, t: str, qidxpos0: int, tidxpos0: int) -> List[Tuple[int, int]]:
    """Filters lonely base pairs that can not stack (both neighbors not in MFE and both neighbors cant pair)"""
    new_bps = []

    for q_index, t_index in bps:
        if not _neighbors_in_mfe(q_index, t_index, bps) and \
           not _neighbors_can_pair(q_index, t_index, q, t, qidxpos0, tidxpos0):
            continue
        new_bps.append((q_index, t_index))
    return new_bps


def filter_lp_mfe(bps: List[Tuple[int, int]], q: str, t: str, qidxpos0: int, tidxpos0: int) -> List[Tuple[int, int]]:
    """Filters lonely base pairs in MFE interaction (both neighbors not in MFE)"""
    new_bps = []

    for q_index, t_index in bps:
        if not _neighbors_in_mfe(q_index, t_index, bps):
            continue
        new_bps.append((q_index, t_index))
    return new_bps


def filter_he(bps: List[Tuple[int, int]], q: str, t: str, qidxpos0: int, tidxpos0: int) -> List[Tuple[int, int]]:
    """Filters helix ends (One neighbor can pair, the other one can not)"""
    new_bps = []

    for q_index, t_index in bps:
        if _neighbors_can_pair(q_index, t_index, q, t, qidxpos0, tidxpos0) == 1:
            continue
        new_bps.append((q_index, t_index))
    return new_bps


def filter_he_mfe(bps: List[Tuple[int, int]], q: str, t: str, qidxpos0: int, tidxpos0: int) -> List[Tuple[int, int]]:
    """Filters helix ends (One neighbor in MFE, the other one is not)"""
    new_bps = []

    for q_index, t_index in bps:
        if _neighbors_in_mfe(q_index, t_index, bps) == 1:
            continue
        new_bps.append((q_index, t_index))
    return new_bps


def get_filter(string: str) -> any:
    candidate_filters = {
        'GU': filter_gu,
        'AU': filter_au,
        'GC': filter_gc,
        'lp': filter_lp,
        'lpMfe': filter_lp_mfe,
        'he': filter_he,
        'heMfe': filter_he_mfe
    }
    return candidate_filters.get(string)

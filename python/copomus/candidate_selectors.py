# Copyright 2019
# Author: Fabio Gutmann <https://github.com/fabio-gut>

from typing import List, Tuple
from copomus.IntaRNA import IntaRNA
from copomus.mutation import parse_mutation_encoding


def get_mfe_bps(query: str, target: str, qidxpos0: int, tidxpos0: int, param_file: str, threads: int, **kwargs) -> \
        List[Tuple[int, int]]:
    """
    Gets a list of all base pairs for a given interaction
    :param query: The query sequence
    :param target: The target sequence
    :param qidxpos0: The counting start index for the query
    :param tidxpos0: The counting start index for the target
    :param param_file: Parameter file used for IntaRNA call
    :param threads: Thread count used by IntaRNA
    :return: List of tuples of indexes in form (query_index, target_index) with each tuple representing a base pair
    """
    candidates = []
    i = IntaRNA()
    data = i.run(query, target, qidxpos0, tidxpos0, 'bpList', threads, param_file=param_file)
    if not data or not data['bpList']:
        return []
    for t in data['bpList'].strip().split(':'):
        t_index, q_index = [int(x) for x in t.strip('(').strip(')').split(',')]
        candidates.append((q_index, t_index))  # note we do query first, then target
    return candidates


def get_mfe_bps_so(query: str, target: str, qidxpos0: int, tidxpos0: int, param_file: str, threads: int, **kwargs) -> \
        List[Tuple[int, int]]:
    """
    Gets a list of all suboptimal base pairs for a given interaction
    :param query: The query sequence
    :param target: The target sequence
    :param qidxpos0: The counting start index for the query
    :param tidxpos0: The counting start index for the target
    :param param_file: Parameter file used for IntaRNA call
    :param threads: Thread count used by IntaRNA
    :return: List of tuples of indexes in form (query_index, target_index) with each tuple representing a base pair
    """
    candidates = []
    i = IntaRNA()
    data = i.run(query, target, qidxpos0, tidxpos0, 'start1,end1,start2,end2', threads, param_file=param_file)
    if not data:
        return []
    data = i.run(query, target, qidxpos0, tidxpos0, 'bpList', threads, 1000, param_file,
                 ['--qRegion', f"{data['start2']}-{data['end2']}", '--tRegion',  f"{data['start1']}-{data['end1']}"])
    seen = set()
    for subopt in data or []:
        for t in subopt['bpList'].strip().split(':'):
            t_index, q_index = [int(x) for x in t.strip('(').strip(')').split(',')]
            if (q_index, t_index) not in seen:
                seen.add((q_index, t_index))
                candidates.append((q_index, t_index))  # note we do query first, then target
    return candidates


def get_mutation_enc(query, target, qidxpos0, tidxpos0, param_file, threads=0, **kwargs):
    """Validate explicit mutations and return their query/target coordinates."""
    return [(q, t) for q, t, _ in parse_mutation_encoding(
        kwargs['enc'], query, target, qidxpos0, tidxpos0)]


def get_selector(string: str) -> any:
    """
    Gets a selector by its string name
    :param string: The name of the selector
    :return: A reference to a selector function
    """
    selectors = {
        'mfe': get_mfe_bps,
        'mfeSO': get_mfe_bps_so,
        'mutEnc': get_mutation_enc
    }
    return selectors.get(string)

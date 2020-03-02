# Copyright 2019
# Author: Fabio Gutmann <https://github.com/fabio-gut>

import sys
from typing import List, Tuple
from copomus.IntaRNA import IntaRNA
from copomus.indexing import idx_to_array_index as idx_aidx


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
        print('No favorable interaction between the specified sequences!')
        sys.exit(1)
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
        print('No favorable interaction between the specified sequences!')
        sys.exit(1)
    data = i.run(query, target, qidxpos0, tidxpos0, 'bpList', threads, 1000, param_file,
                 ['--qRegion', f"{data['start2']}-{data['end2']}", '--tRegion',  f"{data['start1']}-{data['end1']}"])
    for subopt in data:
        for t in subopt['bpList'].strip().split(':'):
            t_index, q_index = [int(x) for x in t.strip('(').strip(')').split(',')]
            if (q_index, t_index) not in candidates:
                candidates.append((q_index, t_index))  # note we do query first, then target
    return candidates


def get_mutation_enc(query: str, target: str, qidxpos0: int, tidxpos0: int, param_file: str, **kwargs) \
        -> List[Tuple[int, int]]:
    """
    Only selects the base pairs specified by the mutation encoding
    :param query: The query sequence
    :param target: The target sequence
    :param qidxpos0: The counting start index for the query
    :param tidxpos0: The counting start index for the target
    :param param_file: Parameter file used for IntaRNA call
    :return: List of tuples of indexes in form (query_index, target_index) with each tuple representing a base pair
    """
    try:
        candidates = []
        enc = kwargs['enc']
        for mut in enc.split(','):  # split block mutations into single mutations
            q, t = mut.split('&')
            q_index = int(q[1:len(q) - 1])  # mutation indices
            t_index = int(t[1:len(t) - 1])

            if query[idx_aidx(q_index, qidxpos0)] != q[0] or target[idx_aidx(t_index, tidxpos0)] != t[0]:
                print('The encoding you specified does not match with the sequences you specified!')
                sys.exit(1)
            candidates.append((q_index, t_index))
        return candidates

    except (IndexError, ValueError, KeyError):
        print('Error parsing the encoding you specified!')
        sys.exit(1)


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

# Copyright 2019
# Author: Fabio Gutmann <https://github.com/fabio-gut>

from copy import deepcopy
from copomus.indexing import idx_to_array_index
from copomus.measures import get_measure
from copomus.mutation import Mutation
from typing import List, Tuple


def gen_flip(query: str, target: str, qidxpos0: int, tidxpos0: int, bp_list: List[Tuple[int, int]], mms: List[str],
             alpha: float, beta: float) -> List[Mutation]:
    mutations = []

    for q_index, t_index in bp_list:
        mm = []
        for m in mms:
            mm.append(get_measure(m, alpha, beta))

        q_index_a = idx_to_array_index(q_index, qidxpos0)
        t_index_a = idx_to_array_index(t_index, tidxpos0)

        bp_mm = f'{target[t_index_a]}{query[q_index_a]}'

        m = Mutation(deepcopy(mm), query, target, qidxpos0, tidxpos0, q_index, t_index, bp_mm)
        mutations.append(m)
    return mutations


def gen_any(query: str, target: str, qidxpos0: int, tidxpos0: int, bp_list: List[Tuple[int, int]], mms: List[str],
            alpha: float, beta: float) -> List[Mutation]:
    mutations = []

    possible_mutations = {
        'AU': ['UA', 'UG', 'CG'],
        'UA': ['AU', 'GU', 'GC'],
        'GC': ['UA', 'UG', 'CG'],
        'CG': ['AU', 'GU', 'GC'],
        'GU': ['UA', 'UG', 'CG'],
        'UG': ['AU', 'GU', 'GC']
    }

    for q_index, t_index in bp_list:
        mm = []
        for m in mms:
            mm.append(get_measure(m, alpha, beta))

        q_index_a = idx_to_array_index(q_index, qidxpos0)
        t_index_a = idx_to_array_index(t_index, tidxpos0)

        pos = possible_mutations[f'{query[q_index_a]}{target[t_index_a]}']

        for p in pos:
            bp_mm = f'{p[0]}{p[1]}'

            m = Mutation(deepcopy(mm), query, target, qidxpos0, tidxpos0, q_index, t_index, bp_mm)
            mutations.append(m)
    return mutations


def gen_specific(query: str, target: str, qidxpos0: int, tidxpos0: int, mms: List[str], alpha: float, beta: float,
                 encoding: str) -> List[Mutation]:
    mutations = []

    mm = []
    for m in mms:
        mm.append(get_measure(m, alpha, beta))

    for mut in encoding.split(','):  # split block mutations into single mutations
        q, t = mut.split('&')
        q_index = int(q[1:len(q)-1])  # mutation indices
        t_index = int(t[1:len(t)-1])

        bp_mm = f'{q[len(q)-1]}{t[len(t)-1]}'

        m = Mutation(deepcopy(mm), query, target, qidxpos0, tidxpos0, q_index, t_index, bp_mm)
        mutations.append(m)
    return mutations


def get_generator(string: str) -> any:
    """
    Gets a mutation generator by its string name
    :param string: The string name
    :return: A reference to a generator function
    """
    generators = {
        'flip': gen_flip,
        'any': gen_any,
        'specific': gen_specific
    }
    return generators.get(string)

# Copyright 2019
# Author: Fabio Gutmann <https://github.com/fabio-gut>

from dataclasses import dataclass, field
from typing import Any
import re
from copomus.indexing import idx_to_array_index


@dataclass()
class Mutation:
    measures: list[Any]
    qw: str  # wildtype sequences
    tw: str
    qidxpos0: int  # indexing start of query
    tidxpos0: int
    query_n: int  # index of mutation on query
    target_n: int
    bp_mm: str  # mm base pair
    bp_ww: str = ''  # ww base pair
    qm: str = ''  # mutated sequences
    tm: str = ''
    ranks: dict = field(default_factory=dict)

    def __post_init__(self):
        mut_pos_q = idx_to_array_index(self.query_n, self.qidxpos0)
        mut_pos_t = idx_to_array_index(self.target_n, self.tidxpos0)
        if not 0 <= mut_pos_q < len(self.qw) or not 0 <= mut_pos_t < len(self.tw):
            raise ValueError('mutation position is outside its sequence')
        if not re.fullmatch('[ACGU]{2}', self.bp_mm):
            raise ValueError('mutated bases must be A, C, G or U')
        self.bp_ww = f'{self.qw[mut_pos_q]}{self.tw[mut_pos_t]}'
        self.qm = self.qw[0:mut_pos_q] + self.bp_mm[0] + self.qw[mut_pos_q+1:]
        self.tm = self.tw[0:mut_pos_t] + self.bp_mm[1] + self.tw[mut_pos_t+1:]

    def __lt__(self, other):
        return self.measures < other.measures

    def __str__(self):
        return f'{self.bp_ww[0]}{self.query_n}{self.bp_mm[0]}&{self.bp_ww[1]}{self.target_n}{self.bp_mm[1]}'


def parse_mutation_encoding(encoding, query, target, qidxpos0, tidxpos0):
    """Validate query-first substitutions without allowing negative array wrapping."""
    mutations = []
    for entry in encoding.split(','):
        match = re.fullmatch(r'([ACGU])(-?\d+)([ACGU])&([ACGU])(-?\d+)([ACGU])', entry)
        if match is None:
            raise ValueError(f'invalid mutation encoding: {entry!r}')
        qw, q, qm, tw, t, tm = match.groups()
        q, t = int(q), int(t)
        qi, ti = idx_to_array_index(q, qidxpos0), idx_to_array_index(t, tidxpos0)
        if not 0 <= qi < len(query) or not 0 <= ti < len(target):
            raise ValueError('mutation position is outside its sequence')
        if query[qi] != qw or target[ti] != tw:
            raise ValueError('mutation encoding does not match the wildtype sequences')
        if qm == qw or tm == tw:
            raise ValueError('both encoded bases must change')
        mutations.append((q, t, qm + tm))
    return mutations

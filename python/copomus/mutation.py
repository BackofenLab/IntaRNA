# Copyright 2019
# Author: Fabio Gutmann <https://github.com/fabio-gut>

from dataclasses import dataclass, field
from typing import List
from copomus.indexing import idx_to_array_index


@dataclass()
class Mutation:
    measures: List[any]
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
        self.bp_ww = f'{self.qw[mut_pos_q]}{self.tw[mut_pos_t]}'
        self.qm = self.qw[0:mut_pos_q] + self.bp_mm[0] + self.qw[mut_pos_q+1:]
        self.tm = self.tw[0:mut_pos_t] + self.bp_mm[1] + self.tw[mut_pos_t+1:]

    def __lt__(self, other):
        return self.measures < other.measures

    def __str__(self):
        return f'{self.bp_ww[0]}{self.query_n}{self.bp_mm[0]}&{self.bp_ww[1]}{self.target_n}{self.bp_mm[1]}'

# Copyright 2019
# Author: Fabio Gutmann <https://github.com/fabio-gut>

import os
from tempfile import gettempdir
from copomus.mutation import Mutation
from copomus.indexing import idx_to_array_index


class Measure:
    def __init__(self):
        self.outcsvcols = []  # list with all needed outcsvcols
        self.out_pipes = []
        self.score = {'ww': None, 'wm': None, 'mw': None, 'mm': None}  # measurement scores

    def calculate_score(self, data: dict, mut_type: str, m: Mutation = None, param_file: str = ''):
        if data:
            self.score[mut_type] = data[self.outcsvcols[0]]

    def post_compute(self):
        pass

    def read_temp_file(self):
        data = []
        for pipe in self.out_pipes:
            file_name = f'CopomuS_{pipe}_{os.getpid()}.temp'
            with open(os.path.join(gettempdir(), file_name), 'r') as f:
                data.append(f.read())
            os.remove(os.path.join(gettempdir(), file_name))
        return data

    # noinspection PyTypeChecker
    def is_valid(self, alpha: float, beta: float = None) -> bool:
        beta = alpha if not beta else beta

        res = (self.score['ww'] < 0 and self.score['mm'] < 0 and
               (self.score['mw'] >= 0 or
                ((self.score['ww'] + alpha < self.score['mw']) and (self.score['mm'] + beta < self.score['mw'])))
               and
               (self.score['wm'] >= 0 or
                ((self.score['ww'] + alpha < self.score['wm']) and (self.score['mm'] + beta < self.score['wm'])))
               )
        return res

    def non_i(self) -> bool:
        """Checks if the Measure as values from non-interactions"""
        return None in self.score.values()

    def __lt__(self, other):
        return False  # overwrite for custom comparisons

    def __eq__(self, other):
        return self.score == other.score


class Energy(Measure):
    def __init__(self, alpha, beta):
        super().__init__()
        self.outcsvcols = ['E']
        self.id = 'E'
        self.alpha = alpha
        self.beta = beta

    def post_compute(self):
        if self.non_i():
            # noinspection PyTypeChecker
            self.score['is_valid'] = 0
            return
        # noinspection PyTypeChecker
        self.score['is_valid'] = int(self.is_valid(self.alpha, self.beta))

    def __lt__(self, other):
        if self.non_i():
            return False
        elif other.non_i():
            return True
        return self.is_valid(self.alpha, self.beta) and not other.is_valid(self.alpha, self.beta)

    def __eq__(self, other):
        if self.non_i() and not other.non_i() or not self.non_i() and other.non_i():
            return False
        if self.non_i() and other.non_i():
            return True
        return self.is_valid(self.alpha, self.beta) == other.is_valid(self.alpha, self.beta)


class MinDeltaEnergy(Measure):
    def __init__(self):
        super().__init__()
        self.outcsvcols = ['E']
        self.id = 'minDeltaE'

    def post_compute(self):
        score = {
            'wm/ww': round(self.score['wm'] - self.score['ww'], 2) if self.score['wm'] and self.score['ww'] else '',
            'mw/ww': round(self.score['mw'] - self.score['ww'], 2) if self.score['mw'] and self.score['ww'] else '',
            'wm/mm': round(self.score['wm'] - self.score['mm'], 2) if self.score['wm'] and self.score['mm'] else '',
            'mw/mm': round(self.score['mw'] - self.score['mm'], 2) if self.score['mw'] and self.score['mm'] else '',
            'ww/mm': round(self.score['ww'] - self.score['mm'], 2) if self.score['ww'] and self.score['mm'] else ''
        }
        self.score = score
        entries = [v for k, v in self.score.items() if v != '' and k != 'ww/mm']
        self.score['min'] = min(entries) if entries else None

    def __lt__(self, other):
        if self.non_i():
            return False
        elif other.non_i():
            return True
        return self.score['min'] > 0 > other.score['min'] or 0 < other.score['min'] < self.score['min']

    def __eq__(self, other):
        if self.non_i() and not other.non_i() or not self.non_i() and other.non_i():
            return False
        if self.non_i() and other.non_i():
            return True
        return self.score['min'] == other.score['min']


class MFECover(Measure):
    def __init__(self):
        super().__init__()
        self.outcsvcols = ['start1', 'end1', 'start2', 'end2']
        self.id = 'mfeCover'

    def post_compute(self):
        if self.non_i():
            # noinspection PyTypeChecker
            self.score['is_valid'] = 0
            return
        # noinspection PyTypeChecker
        self.score['is_valid'] = int(self.in_mfe())

    def in_mfe(self):
        return self.score['ww'] and self.score['mm']

    def calculate_score(self, data: dict, mut_type: str, m: Mutation = None, param_file: str = ''):
        if not data:
            return
        mfe_cover = (m.target_n in range(data['start1'], data['end1']+1)) and \
                    (m.query_n in range(data['start2'], data['end2']+1))
        # noinspection PyTypeChecker
        self.score[mut_type] = mfe_cover

    def __lt__(self, other):
        if self.non_i():
            return False
        elif other.non_i():
            return True
        return self.in_mfe() and not other.in_mfe()

    def __eq__(self, other):
        if self.non_i() and not other.non_i() or not self.non_i() and other.non_i():
            return False
        if self.non_i() and other.non_i():
            return True
        return self.in_mfe() == other.in_mfe()


class EnergyProfileQuery(Measure):
    def __init__(self, alpha, beta):
        super().__init__()
        self.id = 'Eqi'
        self.out_pipes = ['qMinE']
        self.alpha = alpha
        self.beta = beta

    def post_compute(self):
        if self.non_i():
            # noinspection PyTypeChecker
            self.score['is_valid'] = 0
            return
        # noinspection PyTypeChecker
        self.score['is_valid'] = int(self.is_valid(self.alpha, self.beta))

    def calculate_score(self, data: dict, mut_type: str, m: Mutation = None, param_file: str = ''):
        data = self.read_temp_file()[0]
        data = data.strip().split('\n')
        value = data[idx_to_array_index(m.query_n, m.qidxpos0) + 1].split(';')[2]
        if value == 'NA':
            self.score[mut_type] = None
        else:
            # noinspection PyTypeChecker
            self.score[mut_type] = float(value)

    def __lt__(self, other):
        if self.non_i():
            return False
        elif other.non_i():
            return True
        return self.is_valid(self.alpha, self.beta) and not other.is_valid(self.alpha, self.beta)

    def __eq__(self, other):
        if self.non_i() and not other.non_i() or not self.non_i() and other.non_i():
            return False
        if self.non_i() and other.non_i():
            return True
        return self.is_valid(self.alpha, self.beta) == other.is_valid(self.alpha, self.beta)


class EnergyProfileTarget(Measure):
    def __init__(self, alpha, beta):
        super().__init__()
        self.id = 'Eti'
        self.out_pipes = ['tMinE']
        self.alpha = alpha
        self.beta = beta

    def post_compute(self):
        if self.non_i():
            # noinspection PyTypeChecker
            self.score['is_valid'] = 0
            return
        # noinspection PyTypeChecker
        self.score['is_valid'] = int(self.is_valid(self.alpha, self.beta))

    # noinspection PyTypeChecker
    def calculate_score(self, data: dict, mut_type: str, m: Mutation = None, param_file: str = ''):
        data = self.read_temp_file()[0]
        data = data.strip().split('\n')
        value = data[idx_to_array_index(m.target_n, m.tidxpos0) + 1].split(';')[2]
        if value == 'NA':
            self.score[mut_type] = None
        else:
            # noinspection PyTypeChecker
            self.score[mut_type] = float(value)

    def __lt__(self, other):
        if self.non_i():
            return False
        elif other.non_i():
            return True
        return self.is_valid(self.alpha, self.beta) and not other.is_valid(self.alpha, self.beta)

    def __eq__(self, other):
        if self.non_i() and not other.non_i() or not self.non_i() and other.non_i():
            return False
        if self.non_i() and other.non_i():
            return True
        return self.is_valid(self.alpha, self.beta) == other.is_valid(self.alpha, self.beta)


class AccessProfileQuery(Measure):
    def __init__(self):
        super().__init__()
        self.id = 'EDqi'
        self.out_pipes = ['qAcc']

    def post_compute(self):
        if self.non_i():
            # noinspection PyTypeChecker
            self.score['is_valid'] = 0
            return
        # noinspection PyTypeChecker
        self.score['is_valid'] = int(self.is_valid(0))

    def calculate_score(self, data: dict, mut_type: str, m: Mutation = None, param_file: str = ''):
        data = self.read_temp_file()[0]
        data = data.strip().split('\n')
        # noinspection PyTypeChecker
        self.score[mut_type] = float(data[idx_to_array_index(m.query_n, m.qidxpos0) + 2].split('\t')[1])

    # noinspection PyTypeChecker
    def is_valid(self, alpha: float, beta: float = None, gamma: float = 0) -> bool:
        return self.score['ww'] < self.score['mm']

    def __lt__(self, other):
        if self.non_i():
            return False
        elif other.non_i():
            return True
        return (self.is_valid(0) and not other.is_valid(0)) or (
                self.is_valid(0) and other.is_valid(0) and
                self.score['mm'] - self.score['ww'] < other.score['mm'] - other.score['ww'])

    def __eq__(self, other):
        if self.non_i() and not other.non_i() or not self.non_i() and other.non_i():
            return False
        if self.non_i() and other.non_i():
            return True
        return self.is_valid(0) == other.is_valid(0) and \
               self.score['mm'] - self.score['ww'] == other.score['mm'] - other.score['ww']


class ConstantAccessProfileQuery(Measure):
    def __init__(self):
        super().__init__()
        self.id = 'cEDqi'
        self.out_pipes = ['qAcc']

    def post_compute(self):
        if self.non_i():
            # noinspection PyTypeChecker
            self.score['is_valid'] = 0
            return
        # noinspection PyTypeChecker
        self.score['is_valid'] = int(self.is_valid(1))

    def calculate_score(self, data: dict, mut_type: str, m: Mutation = None, param_file: str = ''):
        data = self.read_temp_file()[0]
        data = data.strip().split('\n')
        # noinspection PyTypeChecker
        self.score[mut_type] = float(data[idx_to_array_index(m.query_n, m.qidxpos0) + 2].split('\t')[1])

    # noinspection PyTypeChecker
    def is_valid(self, alpha: float, beta: float = None, gamma: float = 0) -> bool:
        return abs(self.score['ww'] - self.score['mm']) < alpha

    def __lt__(self, other):
        if self.non_i():
            return False
        elif other.non_i():
            return True
        return self.is_valid(1) and not other.is_valid(1)

    def __eq__(self, other):
        if self.non_i() and not other.non_i() or not self.non_i() and other.non_i():
            return False
        if self.non_i() and other.non_i():
            return True
        return self.is_valid(1) == other.is_valid(1)


def get_measure(string: str, alpha: float, beta: float) -> Measure:
    mm = {
        'E': Energy(alpha, beta),
        'minDeltaE': MinDeltaEnergy(),
        'mfeCover': MFECover(),
        'EDqi': AccessProfileQuery(),
        'cEDqi': ConstantAccessProfileQuery(),
        'Eqi': EnergyProfileQuery(alpha, beta),
        'Eti': EnergyProfileTarget(alpha, beta)
    }
    return mm.get(string)

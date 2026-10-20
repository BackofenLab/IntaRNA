#!/usr/bin/env python3
"""Extract numerical gold data from IntaRNA-probabilities tests/examples/*.json.

Source revision: f62778504c8abbec8fc9894da16a357b60707680. Run with the path
of that checkout as argument; redirect stdout to tests/data/seeded-theory.dat.
Values preserve unrounded accessibility; they are NOT native ED gold values.
"""
import json
from decimal import Decimal, getcontext
from pathlib import Path
import sys
getcontext().prec = 50
e = Decimal(1).exp()
def poly(cs):
    return sum(Decimal(c) * e**i for i, c in enumerate(cs))
for path in sorted((Path(sys.argv[1]) / 'tests/examples').glob('*.json')):
    d = json.loads(path.read_text())
    p = d['input']
    n, m = len(p['rna1']), len(p['rna2'])
    def pair(v):
        return f'{v[0]-1} {m-v[1]}' # one-based antiparallel -> internal
    print('CASE', p['id'], n, m, p['max_interaction_length'])
    print('Z', d['ensembles']['seed']['partition_value_at_e'])
    for s in d['seeds']['occurrences']:
        print('S', pair(s['left_pair']), len(s['pairs']))
    denominator = Decimal(d['common_denominator']['value_at_e'])
    for b in d['helpers']['boundaries']:
        print('B', pair(b['left_pair']), pair(b['right_pair']),
              poly(b['accessibility_product_coefficients'])/denominator,
              poly(b['hybrid']['seed']), poly(b['hybrid']['noseed']))
    for a in d['actual_pairs']:
        print('P', pair(a['pair']), int(a['allowed']), a['partition_value_at_e'])
    print('END')

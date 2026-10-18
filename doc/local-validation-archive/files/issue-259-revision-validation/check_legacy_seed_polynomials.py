#!/usr/bin/env python3
"""Independent unit-weight chain-polynomial check of legacy left correction.

Every polynomial coefficient is an integer multiplicity for an explicit pair
chain. This exploratory check is separate from the C++ reverse trace; it does
not prove general correctness or replace its API/oracle regressions.
"""
import json
from collections import defaultdict
from pathlib import Path
from random import Random


def add(out, src, prefix=(), scale=1):
    for chain, weight in src.items():
        out[tuple(sorted(set(chain).union(prefix)))] += weight * scale


def calculate(anchor, seeds, length, no_lp, loop):
    table = {}
    for x in range(anchor[0], -1, -1):
        for y in range(anchor[1], -1, -1):
            start = (x, y)
            current = defaultdict(int)
            if start == anchor:
                current[()] = 1
            elif x < anchor[0] and y < anchor[1]:
                if no_lp:
                    add(current, table[x+1, y+1], [start])
                for k in range(x+1+no_lp, min(anchor[0], x+1+no_lp+loop)+1):
                    for l in range(y+1+no_lp, min(anchor[1], y+1+no_lp+loop)+1):
                        prefix = [start] + ([(x+1, y+1)] if no_lp else [])
                        add(current, table[k, l], prefix)

                def overlaps(a, b):
                    return 0 < b[0]-a[0] <= length-2 and b[0]-a[0] == b[1]-a[1]

                if start in seeds and ((x+length-1 <= anchor[0] and y+length-1 <= anchor[1]) or overlaps(start, anchor)):
                    overlapping = [p for p in seeds if p[0] <= anchor[0] and p[1] <= anchor[1] and overlaps(start, p)]
                    if overlapping:
                        child = min(overlapping)
                        add(current, table[child], [(x+d, y+d) for d in range(child[0]-x)], -1)
                    else:
                        end = (x+length-1, y+length-1)
                        add(current, table[end], [(x+d, y+d) for d in range(length-1)], -1)
                        if no_lp:
                            for k in range(end[0]+1, min(anchor[0], end[0]+loop+1)+1):
                                for l in range(end[1]+1, min(anchor[1], end[1]+loop+1)+1):
                                    add(current, table[k, l], [(x+d, y+d) for d in range(length)], -1)
                    if sum(current.values()) < 0:
                        current = defaultdict(int)
            table[start] = dict(current)
    return table


def check(table):
    for start, polynomial in table.items():
        total = sum(polynomial.values())
        masses = defaultdict(int)
        for chain, coefficient in polynomial.items():
            for pair in chain:
                masses[pair] += coefficient
        if total > 0:
            assert all(0 <= mass <= total for mass in masses.values()), (start, total, dict(masses), polynomial)


def main():
    random = Random(259)
    for trial in range(500):
        anchor = (4, 4)
        length = random.choice([2, 3, 4])
        no_lp = random.choice([False, True])
        loop = random.choice([0, 1, 2, 4])
        seeds = {anchor} | {(i, j) for i in range(5) for j in range(5) if random.random() < 0.3}
        check(calculate(anchor, seeds, length, no_lp, loop))
    result = dict(random_seed=259, trials=500, domain="5x5 all-complementary unit-weight pairs", seed_lengths=[2, 3, 4], loop_limits=[0, 1, 2, 4], no_lp=[False, True], result="All positive-Z cells had pair marginals in [0,Z].", limitation="Exploratory unit-weight polynomial check, not a proof or a C++ trace test.")
    Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()

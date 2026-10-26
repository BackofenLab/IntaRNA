# Copyright 2026
# Author: Alexander Mitrofanov <https://github.com/Alexander-Mitrofanov>

"""Shared CopomuS ranking and tabular output helpers."""

import csv


def calculate_ranks(mutations, measure_ids):
    """Populate total and per-measure ranks after mutations were sorted."""
    rank = 1
    mutations[0].ranks['total'] = rank
    for index in range(1, len(mutations)):
        if mutations[index - 1] < mutations[index]:
            rank += 1
        mutations[index].ranks['total'] = rank

    for measure_index, measure_id in enumerate(measure_ids):
        rank = 1
        measures = [(mutation.measures[measure_index], index)
                    for index, mutation in enumerate(mutations)]
        measures.sort(key=lambda entry: entry[0])
        mutations[measures[0][1]].ranks[measure_id] = rank
        for index in range(1, len(measures)):
            measure, mutation_index = measures[index]
            if measures[index - 1][0] < measure:
                rank += 1
            mutations[mutation_index].ranks[measure_id] = rank


def table(mutations):
    """Return the stable CopomuS header and rows for ranked mutations."""
    headers = [column for measure in mutations[0].measures
               for column in ([f'{measure.id}_rank']
                              + [f'{measure.id}_{key}' for key in measure.score])]
    header = ['mutation', 'rank', 'qIndex', 'tIndex', 'bpWildtype', 'bpMutated', *headers]
    rows = []
    for mutation in mutations:
        values = []
        for measure in mutation.measures:
            values.append(mutation.ranks[measure.id])
            values.extend('NA' if value is None else value for value in measure.score.values())
        rows.append([str(mutation), mutation.ranks['total'], mutation.query_n,
                     mutation.target_n, mutation.bp_ww, mutation.bp_mm, *values])
    return header, rows


def write_csv(mutations, handle, delimiter='\t'):
    """Write ranked mutations with deterministic UTF-8/LF-compatible CSV."""
    writer = csv.writer(handle, delimiter=delimiter, lineterminator='\n')
    header, rows = table(mutations)
    writer.writerow(header)
    writer.writerows(rows)

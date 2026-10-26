# Copyright 2026
# Author: Alexander Mitrofanov <https://github.com/Alexander-Mitrofanov>

"""Asynchronous CopomuS workflow for a JSON-to-WASM browser boundary."""

import io
import json
import math

from copomus.IntaRNA import IntaRNA
from copomus.candidate_filters import get_filter
from copomus.mutation import parse_mutation_encoding
from copomus.mutation_generators import get_generator
from copomus.workflow import calculate_ranks, table, write_csv


MEASURES = ('E', 'minDeltaE', 'mfeCover', 'EDqi', 'cEDqi', 'Eqi', 'Eti')
SELECTORS = ('mfe', 'mfeSO')
FILTERS = ('GU', 'AU', 'GC', 'lp', 'lpMfe', 'he', 'heMfe')
GENERATORS = ('flip', 'any')
MAX_CANDIDATES = 100
CONTEXT_COLUMNS = ('id1', 'start1', 'end1', 'id2', 'start2', 'end2',
                   'subseqDP', 'hybridDP', 'E')


async def _call(predict, query, target, request, *, columns, number=1,
                outputs=(), parameters=None):
    payload = {
        'query': query,
        'target': target,
        'qIdxPos0': request['qIdxPos0'],
        'tIdxPos0': request['tIdxPos0'],
        'columns': list(dict.fromkeys(columns)),
        'outNumber': number,
        'outputs': list(dict.fromkeys(outputs)),
        'parameters': parameters or {},
    }
    response = json.loads(await predict(json.dumps(payload, separators=(',', ':'))))
    if response.get('error'):
        raise RuntimeError(response['error'])
    rows = IntaRNA.output_to_dict(response.get('csv', ''))
    return rows, response.get('files', {})


def _base_pairs(rows):
    candidates = []
    seen = set()
    for row in rows or []:
        for pair in (row.get('bpList') or '').split(':'):
            if not pair:
                continue
            target_index, query_index = [int(value) for value in pair.strip('()').split(',')]
            candidate = (query_index, target_index)
            if candidate not in seen:
                seen.add(candidate)
                candidates.append(candidate)
    return candidates


async def _select_candidates(request, predict):
    options = request['copomus']
    if options['mutationEncoding']:
        return [(query, target) for query, target, _ in parse_mutation_encoding(
            options['mutationEncoding'], request['query'], request['target'],
            request['qIdxPos0'], request['tIdxPos0'])]

    if options['candidateSelection'] == 'mfe':
        rows, _ = await _call(predict, request['query'], request['target'], request,
                              columns=['bpList'])
        return _base_pairs(rows)

    optimum, _ = await _call(
        predict, request['query'], request['target'], request,
        columns=['start1', 'end1', 'start2', 'end2'])
    if not optimum:
        return []
    optimum = optimum[0]
    rows, _ = await _call(
        predict, request['query'], request['target'], request,
        columns=['bpList'], number=1000,
        parameters={'qRegion': f"{optimum['start2']}-{optimum['end2']}",
                    'tRegion': f"{optimum['start1']}-{optimum['end1']}"})
    return _base_pairs(rows)


def _validate(request):
    required = {'query', 'target', 'qIdxPos0', 'tIdxPos0', 'copomus'}
    if not isinstance(request, dict) or not required <= request.keys():
        raise ValueError('Incomplete mutation-design request')
    options = request['copomus']
    if (not isinstance(options, dict)
            or not isinstance(options.get('measures'), list)
            or not options['measures']
            or len(set(options['measures'])) != len(options['measures'])
            or any(measure not in MEASURES for measure in options['measures'])):
        raise ValueError('Choose one or more available CopomuS measures')
    if options.get('candidateSelection') not in SELECTORS:
        raise ValueError('Choose an available candidate-selection method')
    if (not isinstance(options.get('candidateFilters'), list)
            or len(set(options['candidateFilters'])) != len(options['candidateFilters'])
            or any(candidate_filter not in FILTERS
                   for candidate_filter in options['candidateFilters'])):
        raise ValueError('Choose available candidate filters without duplicates')
    if options.get('generator') not in GENERATORS:
        raise ValueError('Choose an available mutation generator')
    if not isinstance(options.get('mutationEncoding'), str):
        raise ValueError('Mutation encoding must be text')
    if not all(isinstance(options.get(key), (int, float))
               and math.isfinite(options[key]) for key in ('alpha', 'beta')):
        raise ValueError('alpha and beta must be finite numbers')


async def run(request_json, predict, progress=lambda _message: None):
    """Run CopomuS using an awaitable JSON predictor supplied by JavaScript."""
    request = json.loads(request_json)
    _validate(request)
    options = request['copomus']
    progress('Selecting compensatory-mutation candidates…')
    base_pairs = await _select_candidates(request, predict)
    if not options['mutationEncoding']:
        for candidate_filter in options['candidateFilters']:
            base_pairs = get_filter(candidate_filter)(
                base_pairs, request['query'], request['target'],
                request['qIdxPos0'], request['tIdxPos0'])
    if not base_pairs:
        raise ValueError('There are no base-pair candidates left after filtering')

    if options['mutationEncoding']:
        mutations = get_generator('specific')(
            request['query'], request['target'], request['qIdxPos0'],
            request['tIdxPos0'], options['measures'], options['alpha'],
            options['beta'], options['mutationEncoding'])
    else:
        mutations = get_generator(options['generator'])(
            request['query'], request['target'], request['qIdxPos0'],
            request['tIdxPos0'], base_pairs, options['measures'],
            options['alpha'], options['beta'])
    if len(mutations) > MAX_CANDIDATES:
        raise ValueError(
            f'This browser run produced {len(mutations)} mutations; limit the '
            f'candidate set to at most {MAX_CANDIDATES} with filters or an explicit encoding')

    wildtype = None
    for mutation_index, mutation in enumerate(mutations, 1):
        progress(f'Evaluating mutation {mutation_index} of {len(mutations)}…')
        columns = list(dict.fromkeys(
            [column for measure in mutation.measures for column in measure.outcsvcols]
            + list(CONTEXT_COLUMNS)))
        outputs = list(dict.fromkeys(
            output for measure in mutation.measures for output in measure.out_pipes))
        mutation_visualizations = {}
        for mode in ('ww', 'wm', 'mw', 'mm'):
            if mode == 'ww' and wildtype is not None:
                rows, profiles = wildtype
            else:
                query = mutation.qw if mode[0] == 'w' else mutation.qm
                target = mutation.tw if mode[1] == 'w' else mutation.tm
                rows, profiles = await _call(
                    predict, query, target, request, columns=columns, outputs=outputs)
                if mode == 'ww':
                    wildtype = rows, profiles
            data = rows[0] if rows else False
            mutation_visualizations[mode] = data or None
            for measure in mutation.measures:
                measure.profiles = profiles
                measure.calculate_score(data, mode, mutation, '')
        for measure in mutation.measures:
            measure.post_compute()
        mutation.browser_visualizations = mutation_visualizations

    mutations.sort()
    calculate_ranks(mutations, options['measures'])
    header, rows = table(mutations)
    output = io.StringIO()
    write_csv(mutations, output)
    candidates = []
    for mutation, values in zip(mutations, rows):
        candidates.append({
            'values': dict(zip(header, values)),
            'queryWildtype': mutation.qw,
            'targetWildtype': mutation.tw,
            'queryMutated': mutation.qm,
            'targetMutated': mutation.tm,
            'interactions': mutation.browser_visualizations,
        })
    return json.dumps({'tsv': output.getvalue(), 'columns': header,
                       'candidates': candidates}, separators=(',', ':'))

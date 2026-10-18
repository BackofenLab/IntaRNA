#!/usr/bin/env python3
"""Audit existing IntaRNA profiling evidence; never runs a benchmark.

Run after benchmark.py, supplement.py, and structural_checks.py finish:
    python3 profiling/scripts/audit.py
The default expectation is 47 original + 9 supplementary cases, each with one
warmup and three measurements. --allow-incomplete permits partial coverage, but
does not suppress evidence-integrity failures. All findings are preserved in
results/audit.json and results/audit.md. Python standard library only.
"""
import argparse
import csv
import hashlib
import io
import json
import math
import statistics
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from benchmark import canonical

SUPPLEMENTARY = {
    'cache_computed_control', 'cache_reused', 'blas_one_control', 'blas_unset',
    'blas_twelve', 'memory_s', 'memory_s_window300',
    'cache_matched_computed', 'cache_matched_reused',
}
HEADER = ['id1', 'id2', 'start1', 'end1', 'start2', 'end2', 'E']
CORRECTED_FIELDS = {'csv_valid', 'canonical_sha256', 'output_rows'}
TIME_KEYS = ['wall_s', 'user_s', 'system_s', 'max_rss_kib', 'major_faults',
             'minor_faults', 'voluntary_context_switches',
             'involuntary_context_switches', 'fs_inputs', 'fs_outputs',
             'time_exit_code']


def digest(data):
    return hashlib.sha256(data).hexdigest()


def identity(row):
    return (row.get('case'), row.get('phase'), row.get('repetition'))


def same(a, b):
    if isinstance(a, (int, float)) and not isinstance(a, bool):
        return isinstance(b, (int, float)) and math.isclose(a, b, rel_tol=1e-11, abs_tol=1e-11)
    return a == b


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--expected-original', type=int, default=47)
    parser.add_argument('--expected-supplementary', type=int, default=9)
    parser.add_argument('--repetitions', type=int, default=3)
    parser.add_argument('--allow-incomplete', action='store_true')
    args = parser.parse_args()
    root = args.root.resolve()
    result = root / 'results'
    errors, warnings, incompleteness = [], [], []

    def problem(message, incomplete=False):
        (incompleteness if incomplete else errors).append(message)

    def read_json(path, default):
        try:
            return json.loads(path.read_text())
        except (OSError, ValueError) as exc:
            problem(f'{path.relative_to(root)}: {exc}', incomplete=not path.exists())
            return default

    def read_jsonl(path):
        try:
            lines = path.read_text().splitlines()
        except OSError as exc:
            problem(f'{path.relative_to(root)}: {exc}', incomplete=True)
            return []
        rows = []
        for number, line in enumerate(lines, 1):
            try:
                rows.append(json.loads(line))
            except ValueError as exc:
                problem(f'{path.name}:{number}: invalid JSON: {exc}')
        return rows

    cases = read_json(root / 'cases.json', [])
    definitions = {case['name']: case for case in cases}
    if len(definitions) != len(cases):
        problem('Duplicate names in cases.json')
    original = [name for name in definitions if name not in SUPPLEMENTARY]
    supplementary = [name for name in definitions if name in SUPPLEMENTARY]
    for label, observed, expected in [('original', len(original), args.expected_original),
                                       ('supplementary', len(supplementary), args.expected_supplementary)]:
        if observed != expected:
            problem(f'Case definition coverage: {observed} {label}, expected {expected}', incomplete=True)

    raw = read_jsonl(result / 'runs.jsonl')
    validated = read_jsonl(result / 'validated-runs.jsonl')
    summaries = read_json(result / 'summary.json', [])
    failures_file = read_json(result / 'failures.json', [])
    environment = read_json(root / 'environment.json', {})
    expected_binary = environment.get('binary_sha256')
    binary_path = root / 'source/src/bin/IntaRNA'
    if binary_path.exists():
        if digest(binary_path.read_bytes()) != expected_binary:
            problem('Current executable SHA-256 differs from environment.json')
    else:
        warnings.append('Current executable unavailable; saved per-run provenance still audited')

    manifest_checks = []
    for manifest_path in [root / 'inputs/manifest.json', root / 'inputs/cache/manifest.json']:
        manifest = read_json(manifest_path, {})
        for name, metadata in manifest.items():
            path = manifest_path.parent / name
            valid = path.exists() and digest(path.read_bytes()) == metadata.get('sha256')
            manifest_checks.append({'path': str(path.relative_to(root)), 'hash_valid': valid})
            if not valid:
                problem(f'Input manifest SHA-256 mismatch/missing file: {path.relative_to(root)}')

    raw_keys = Counter(identity(row) for row in raw)
    validated_keys = Counter(identity(row) for row in validated)
    for key, count in raw_keys.items():
        if count != 1:
            problem(f'Duplicate raw record {key}: {count}')
    if raw_keys != validated_keys:
        problem('validated-runs.jsonl does not contain exactly the raw run identities')
    validated_by_key = {identity(row): row for row in validated}
    corrected = []
    corrections = []
    stderr_nonempty = []
    timing_discrepancies = []
    for record in raw:
        key = identity(record)
        name, phase, repetition = key
        label = f'{name}.{phase}{repetition}'
        if name not in definitions:
            problem(f'Unknown case in raw evidence: {name}')
        elif record.get('command', [])[-len(definitions[name]['args']):] != definitions[name]['args']:
            problem(f'{label}: recorded command differs from case definition')
        if record.get('env_overrides', {}) != definitions.get(name, {}).get('env_overrides', {}):
            problem(f'{label}: environment override differs from case definition')
        if record.get('binary_sha256') != expected_binary:
            problem(f'{label}: binary SHA-256 differs from saved environment')
        if not isinstance(record.get('elapsed_s'), (int, float)) or record.get('elapsed_s', 0) <= 0:
            problem(f'{label}: invalid elapsed time')
        try:
            started = datetime.fromisoformat(record['started_utc'])
            ended = datetime.fromisoformat(record['ended_utc'])
            if ended < started:
                problem(f'{label}: end timestamp precedes start')
        except (KeyError, ValueError, TypeError):
            problem(f'{label}: missing/invalid UTC timestamps')

        stdout_path = result / (label + '.stdout')
        try:
            output = stdout_path.read_bytes()
        except OSError as exc:
            problem(f'{label}: missing stdout: {exc}')
            continue
        if digest(output) != record.get('output_sha256'):
            problem(f'{label}: raw stdout SHA-256 mismatch')
        normalized = canonical(output)
        try:
            rows = list(csv.reader(io.StringIO(normalized.decode()), delimiter=';'))
            csv_valid = bool(rows) and rows[0] == HEADER and all(len(row) == len(HEADER) for row in rows)
            numeric_valid = csv_valid and all(
                all(int(value) != 0 for value in row[2:6]) and math.isfinite(float(row[6]))
                for row in rows[1:])
        except (ValueError, UnicodeError):
            rows, csv_valid, numeric_valid = [], False, False
        derived = dict(record, canonical_sha256=digest(normalized), csv_valid=csv_valid,
                       output_rows=max(0, len(rows) - 1))
        corrected.append(derived)
        delta = {field: {'raw': record.get(field), 'derived': derived[field]}
                 for field in CORRECTED_FIELDS if record.get(field) != derived[field]}
        if delta:
            corrections.append({'case': name, 'phase': phase, 'repetition': repetition,
                                'fields': delta, 'comment_lines': sum(line.startswith(b'#') for line in output.splitlines())})
        if record.get('exit_code') == 0 and not numeric_valid:
            problem(f'{label}: successful command produced invalid CSV schema or numeric fields')
        if key in validated_by_key and validated_by_key[key] != derived:
            problem(f'{label}: validated record differs from independent derivation (or altered immutable fields)')

        stderr_path = result / (label + '.stderr')
        if not stderr_path.exists():
            problem(f'{label}: missing stderr file')
        elif stderr_path.stat().st_size:
            stderr_nonempty.append({'run': label, 'bytes': stderr_path.stat().st_size,
                                    'preview': stderr_path.read_text(errors='replace')[:600]})
        time_path = result / (label + '.time')
        if not time_path.exists():
            problem(f'{label}: missing GNU time file')
            continue
        values = time_path.read_text().strip().splitlines()
        if values and len(values[-1].split('\t')) == len(TIME_KEYS):
            try:
                time_fields = dict(zip(TIME_KEYS, map(float, values[-1].split('\t'))))
                for field, value in time_fields.items():
                    if not same(value, record.get(field)):
                        problem(f'{label}: {field} differs from GNU time evidence')
                if record.get('exit_code') == 0 and time_fields['time_exit_code'] != 0:
                    problem(f'{label}: GNU time and harness exit codes disagree')
                difference = record['elapsed_s'] - time_fields['wall_s']
                timing_discrepancies.append({'run': label, 'python_minus_gnu_s': difference})
                if difference < -0.02:
                    problem(f'{label}: monotonic wall time implausibly below GNU time')
            except (ValueError, TypeError):
                problem(f'{label}: invalid GNU time data')
        elif record.get('exit_code') == 0:
            problem(f'{label}: successful run lacks parseable GNU time resources')

    by_case = defaultdict(list)
    for row in corrected:
        by_case[row['case']].append(row)
    expected_slots = {('warmup', 0)} | {('measure', i) for i in range(1, args.repetitions + 1)}
    coverage = []
    scatter = []
    recomputed = {}
    for name, definition in definitions.items():
        records = by_case[name]
        observed = Counter((row['phase'], row['repetition']) for row in records)
        missing = sorted(expected_slots - observed.keys())
        extra = sorted(observed.keys() - expected_slots)
        successful = [row for row in records if row['exit_code'] == 0]
        measured = [row for row in successful if row['phase'] == 'measure']
        status = {'case': name, 'warmups': sum(row['phase'] == 'warmup' for row in records),
                  'measurements': sum(row['phase'] == 'measure' for row in records),
                  'successful_measurements': len(measured), 'missing': missing, 'extra': extra,
                  'failed_runs': [list(identity(row)) for row in records if row['exit_code'] != 0],
                  'timeouts': sum(bool(row.get('timeout')) for row in records),
                  'stable_all_successful_outputs': len({row['canonical_sha256'] for row in successful}) <= 1}
        coverage.append(status)
        if missing:
            problem(f'{name}: missing slots {missing}', incomplete=True)
        if extra:
            problem(f'{name}: unexpected slots {extra}')
        if not status['stable_all_successful_outputs']:
            problem(f'{name}: successful warmup/measurement outputs differ')
        if not measured:
            continue
        times = [row['elapsed_s'] for row in measured]
        med = statistics.median(times)
        entry = {'case': name, 'category': definition['category'], 'n': len(measured),
                 'median_s': med, 'min_s': min(times), 'max_s': max(times),
                 'stable_output': len({row['canonical_sha256'] for row in measured}) == 1,
                 'rows': measured[0]['output_rows'], 'canonical_sha256': measured[0]['canonical_sha256']}
        for source, target, divisor, aggregate in [
            ('user_s', 'median_user_s', 1, statistics.median),
            ('system_s', 'median_system_s', 1, statistics.median),
            ('max_rss_kib', 'median_rss_mib', 1024, statistics.median),
            ('max_rss_kib', 'peak_rss_mib', 1024, max)]:
            if all(source in row for row in measured):
                entry[target] = aggregate(row[source] / divisor for row in measured)
        recomputed[name] = entry
        scatter.append({'case': name, 'n': len(times), 'times_s': times, 'median_s': med,
                        'range_s': max(times) - min(times),
                        'range_over_median_pct': 100 * (max(times) - min(times)) / med,
                        'sample_cv_pct': 100 * statistics.stdev(times) / statistics.mean(times) if len(times) > 1 else None,
                        'max_deviation_from_median_pct': 100 * max(abs(t - med) for t in times) / med,
                        'median_cpu_per_wall': statistics.median((r.get('user_s', 0) + r.get('system_s', 0)) / r['elapsed_s'] for r in measured)})

    saved_summary = {row['case']: row for row in summaries}
    if len(saved_summary) != len(summaries):
        problem('Duplicate cases in summary.json')
    if saved_summary.keys() != recomputed.keys():
        problem('summary.json case coverage differs from successful measured runs')
    for name, entry in recomputed.items():
        for field, value in entry.items():
            if not same(value, saved_summary.get(name, {}).get(field)):
                problem(f'{name}: summary field {field} differs from independent recomputation')
    csv_path = result / 'summary.csv'
    if csv_path.exists():
        csv_summary = list(csv.DictReader(csv_path.open()))
        if len(csv_summary) != len(summaries):
            problem('summary.csv and summary.json row counts differ')
        for row in csv_summary:
            reference = saved_summary.get(row.get('case'), {})
            for field, value in reference.items():
                if row.get(field) != str(value):
                    problem(f'{row.get("case")}: summary.csv field {field} differs from JSON summary')
    else:
        problem('summary.csv missing', incomplete=True)
    failed = [row for row in corrected if row['exit_code'] != 0 or not row['csv_valid']]
    if Counter(identity(row) for row in failed) != Counter(identity(row) for row in failures_file):
        problem('failures.json differs from independently derived failures')
    if failed:
        warnings.append(f'{len(failed)} command/schema failures retained; successful-only summaries require coverage caveat')

    groups = []
    def equality(label, names, required=False, performance_comparison_valid=True, caveat=None):
        present = {name: sorted({r['canonical_sha256'] for r in by_case[name]
                                if r['exit_code'] == 0 and r['csv_valid']}) for name in names if name in definitions}
        complete = len(present) == len(names) and all(present.values())
        equal = complete and len({value for hashes in present.values() for value in hashes}) == 1
        groups.append({'group': label, 'cases': names, 'hashes': present,
                       'complete': complete, 'equal': equal, 'required': required,
                       'performance_comparison_valid': performance_comparison_valid,
                       'caveat': caveat})
        if complete and not equal:
            (errors if required else warnings).append(f'{label}: canonical outputs differ')
    equality('thread_scaling', ['throughput_1t', 'throughput_2t', 'throughput_4t', 'throughput_6t', 'throughput_12t'], True)
    equality('prediction_windows', ['length_t3000_q100', 'window_300', 'window_600'])
    equality('memory_windows_model_x', ['memory_t3000_q800', 'memory_window300'])
    equality('memory_windows_model_s', ['memory_s', 'memory_s_window300'])
    equality('cached_accessibility_original_observational',
             ['accessibility_computed', 'cache_computed_control', 'cache_reused'],
             performance_comparison_valid=False,
             caveat='ED reload reduces maximum available length by one for dangling-end treatment '
                    '(AccessibilityFromStream.cpp:90). Original computed/reused effective caps differ '
                    '(target/query 150/100 versus 149/99). Equal observed outputs do not remove this confound.')
    equality('cached_accessibility_matched_caps', ['cache_matched_computed', 'cache_matched_reused'], True,
             caveat='Both commands explicitly cap target/query interaction lengths at 149/99; '
                    'cache generation remains a separate one-time cost.')
    equality('blas_environment', ['bio_fhla_oxys', 'blas_one_control', 'blas_unset', 'blas_twelve'], True)
    equality('duplicate_default_case', ['length_t100_q100', 'algorithm_default'], True)

    matched_cases = [definitions.get(name) for name in ['cache_matched_computed', 'cache_matched_reused']]
    for case in matched_cases:
        if case is not None:
            if any(arg.startswith('--intLenMax=') for arg in case['args']) or not {'--qIntLenMax=99', '--tIntLenMax=149'}.issubset(case['args']):
                problem(f"{case['name']}: matched cache comparison lacks explicit common 149/99 caps")

    structural = read_json(result / 'structural/checks.json', {})
    structural_hashes = {}
    for row in structural.get('runs', []):
        path = result / 'structural' / (row['case'] + '.stdout')
        if not path.exists():
            problem(f"Structural output missing for {row['case']}")
            continue
        data = canonical(path.read_bytes())
        structural_hashes[row['case']] = digest(data)
        if digest(data) != row.get('canonical_sha256'):
            problem(f"Structural stdout hash mismatch for {row['case']}")
        if row.get('exit_code') != 0:
            problem(f"Structural command failed for {row['case']}")
        cells = list(csv.reader(io.StringIO(data.decode()), delimiter=';'))
        if not cells or cells[0] != HEADER + ['bpList'] or any(len(cell) != 8 for cell in cells):
            problem(f"Structural CSV invalid for {row['case']}")
    structural_comparisons = {}
    for label, first, second, required in [
        ('threads_1_vs_12', 'throughput_1t', 'throughput_12t', True),
        ('cached_vs_computed', 'cache_computed_control', 'cache_reused', False),
        ('cache_matched_caps', 'cache_matched_computed', 'cache_matched_reused', True)]:
        equal = first in structural_hashes and second in structural_hashes and structural_hashes[first] == structural_hashes[second]
        structural_comparisons[label] = {'equal': equal, 'required': required,
                                        'performance_comparison_valid': label != 'cached_vs_computed'}
        if structural and structural.get('comparisons', {}).get(label) != equal:
            problem(f'Structural comparison summary mismatch: {label}')
        if required and first in structural_hashes and second in structural_hashes and not equal:
            problem(f'Structural base-pair outputs differ: {label}')

    thread_scaling = []
    baseline = recomputed.get('throughput_1t', {}).get('median_s')
    if baseline:
        for threads in [1, 2, 4, 6, 12]:
            item = recomputed.get(f'throughput_{threads}t')
            if item:
                speedup = baseline / item['median_s']
                thread_scaling.append({'threads': threads, 'median_s': item['median_s'],
                                       'speedup': speedup, 'parallel_efficiency': speedup / threads,
                                       'predictions_per_s': 48 / item['median_s']})

    errors = list(dict.fromkeys(errors))
    warnings = list(dict.fromkeys(warnings))
    incompleteness = list(dict.fromkeys(incompleteness))
    report = {'audited_at_utc': datetime.now(timezone.utc).isoformat(),
              'status': 'fail' if errors or (incompleteness and not args.allow_incomplete) else ('partial' if incompleteness else 'pass'),
              'expected_original_cases': args.expected_original,
              'expected_supplementary_cases': args.expected_supplementary,
              'defined_cases': len(cases), 'raw_records': len(raw), 'validated_records': len(validated),
              'expected_records': (args.expected_original + args.expected_supplementary) * (args.repetitions + 1),
              'errors': errors, 'warnings': warnings, 'incompleteness': incompleteness,
              'input_hash_checks': manifest_checks, 'coverage': coverage,
              'raw_metadata_corrections': corrections, 'failures': failed,
              'nonempty_stderr': stderr_nonempty, 'output_equality_groups': groups,
              'structural_comparisons': structural_comparisons,
              'scatter': sorted(scatter, key=lambda x: x['range_over_median_pct'], reverse=True),
              'thread_scaling': thread_scaling,
              'largest_wrapper_time_differences': sorted(timing_discrepancies, key=lambda x: abs(x['python_minus_gnu_s']), reverse=True)[:10],
              'notes': ['Canonicalization uses the current benchmark.canonical and ignores # comment lines.',
                        'Raw records and benchmark evidence are never edited by this audit.',
                        'Window equality is observed, not presumed; cache/BLAS/thread equality is checked explicitly.',
                        'Only the matched 149/99 cache pair supports a cache performance comparison; original cache timings are confounded.',
                        'Three samples support descriptive scatter, not tight confidence intervals.',
                        'Throughput uses 48 attempted target-query predictions, not reported CSV row count.']}
    result.mkdir(exist_ok=True)
    (result / 'audit.json').write_text(json.dumps(report, indent=2) + '\n')
    lines = ['# Profiling evidence audit', '', f"Status: **{report['status'].upper()}**", '',
             f"{len(cases)} case definitions; {len(raw)}/{report['expected_records']} raw records; "
             f"{len(validated)} validated records; {len(failed)} command/schema failures.", '',
             f'{len(corrections)} raw collection-time metadata records corrected from preserved stdout. '
             'Original raw records remain unchanged.', '', '## Findings', '']
    lines += [f'- ERROR: {message}' for message in errors]
    lines += [f'- INCOMPLETE: {message}' for message in incompleteness]
    lines += [f'- Note: {message}' for message in warnings]
    if not errors and not incompleteness and not warnings:
        lines.append('- No integrity, coverage, or equality issues found.')
    lines += ['', '## Output comparisons', '', '| Comparison | Complete | Equal | Required | Performance comparison valid |', '|---|---|---|---|---|']
    lines += [f"| {g['group']} | {g['complete']} | {g['equal']} | {g['required']} | {g['performance_comparison_valid']} |" for g in groups]
    lines += ['', 'Original cached/computed timings are observational: their effective maximum interaction lengths differ by one. '
              'Only the explicit target/query 149/99 matched pair supports the cache performance comparison. '
              'The JSON companion separately checks base-pair-list output equality.', '']
    lines += ['', '## Largest relative timing ranges', '', '| Case | n | Median s | Min–max s | Range / median |', '|---|---:|---:|---:|---:|']
    for item in report['scatter'][:15]:
        lines.append(f"| {item['case']} | {item['n']} | {item['median_s']:.6f} | "
                     f"{min(item['times_s']):.6f}–{max(item['times_s']):.6f} | {item['range_over_median_pct']:.1f}% |")
    lines += ['', '## Coverage', '', '| Case | Warmup | Measured | Successful measured | Timeouts | Stable outputs |', '|---|---:|---:|---:|---:|---|']
    lines += [f"| {c['case']} | {c['warmups']} | {c['measurements']} | {c['successful_measurements']} | {c['timeouts']} | {c['stable_all_successful_outputs']} |" for c in coverage]
    lines += ['', 'The JSON companion contains detailed provenance, corrections, scatter, and file checks.', '']
    (result / 'audit.md').write_text('\n'.join(lines))
    print(f"Audit {report['status']}: {len(cases)} cases, {len(raw)} records, {len(errors)} integrity errors, "
          f'{len(incompleteness)} incomplete-coverage findings, {len(warnings)} warnings')
    print(result / 'audit.md')
    return 1 if report['status'] == 'fail' else 0


if __name__ == '__main__':
    sys.exit(main())

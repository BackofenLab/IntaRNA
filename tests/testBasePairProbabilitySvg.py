"""SVG semantics, original coordinates, seed provenance and output lifecycle."""
import csv
import gzip
import math
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

BIN = sys.argv[1]
COMMON = ['--energy=B', '--acc=N', '--model=P', '--mode=M', '--seedBP=2',
          '--threads=1', '--outNoLP=false', '--outNoGUend=false', '--outNumber=0',
          '--default-log-file=/dev/null']
NS = {'s': 'http://www.w3.org/2000/svg'}


def run(t='GGAGGG', q='CCCC', options=(), success=True):
    keys = {o.split('=')[0] for o in options}
    if keys & {'--tAcc', '--qAcc'}:
        keys.add('--acc')
    p = subprocess.run([BIN, f'--target={t}', f'--query={q}',
                        *[o for o in COMMON if o.split('=')[0] not in keys], *options],
                       text=True, capture_output=True)
    assert (p.returncode == 0) == success, (options, p.stdout, p.stderr)
    return p


def load(path):
    return ET.fromstring(gzip.decompress(path.read_bytes()) if str(path).endswith('.gz') else path.read_bytes())


def pairs(root):
    return {(int(r.attrib['data-target']), int(r.attrib['data-query'])): r
            for r in root.findall('.//s:rect[@data-type="base-pair"]', NS)}


def seeds(root):
    return {(int(r.attrib['data-target']), int(r.attrib['data-query']))
            for r in root.findall('.//s:rect[@class="seed"][@data-target]', NS)}


with tempfile.TemporaryDirectory(prefix='intarna-bpsvg-') as tmp:
    d = Path(tmp)
    svg, table = d/'pairs.svg', d/'pairs.csv'

    def plot(t='GGAGGG', q='CCCC', opts=()):
        run(t, q, [f'--out=bpsvg:{svg}', *opts])
        return load(svg)

    root = plot(opts=[f'--out=bpProb:{table}'])
    cells = pairs(root)
    rows = list(csv.reader(table.open(), delimiter=';'))
    assert len(cells) == 24
    assert '+1' in [e.text for e in root.findall('.//s:text[@class="index"]', NS)]
    assert not root.findall('.//s:line[@data-index="1"]', NS)
    assert root.find('.//s:text[@class="subtitle"]', NS) is None
    query_axis = root.find('.//s:text[@class="axis"][@x]', NS)
    assert query_axis.text == "query (5' to 3')"
    assert float(query_axis.attrib['y']) > max(float(r.attrib['y']) for r in cells.values())+20
    # Original target coordinate 1 is at the bottom; query coordinate 1 at left.
    for i in range(1, 7):
        for j in range(1, 5):
            assert float(cells[i, j].attrib['y']) == float(cells[1, 1].attrib['y'])-20*(i-1)
            assert float(cells[i, j].attrib['x']) == float(cells[1, 1].attrib['x'])+20*(j-1)
    for i, row in enumerate(rows[1:], 1):
        for j, value in enumerate(row[1:], 1):
            r = cells[i, j]
            assert float(r.attrib['data-probability']) == float(value)
            bounds = [0.01, 0.1, 0.25, 0.5, 0.75, 0.9]
            expected_color = sum(float(value) >= boundary for boundary in bounds)
            assert f'p{expected_color}' in r.attrib['class'].split()
            title = r.find('s:title', NS).text
            assert title.startswith(f'Base-pair probability (target {i}, query {j}): ')
            formatted = f'{float(value):.2e}' if float(value) < .001 else f'{float(value):.3f}'
            assert title.endswith(': '+formatted), title
    assert all(float(cells[3, j].attrib['data-probability']) == 0 for j in range(1, 5))
    frame = root.findall('.//s:rect[@data-type="intramolecular-pairing"]', NS)
    assert len(frame) == 2*(6+4)
    assert all(float(r.attrib['data-probability']) == 0 for r in frame)
    assert all('Intramolecular pairing probability' in r.find('s:title', NS).text for r in frame)
    css = root.find('s:style', NS).text
    assert all(f'.p{i} {{ fill:' in css for i in range(7))
    assert len(root.findall('.//s:g[@class="probability-legend"]/s:rect', NS)) == 7
    # All admitted computed two-pair stacks, in original query orientation.
    expected = {(i+k+1, 4-j-k) for i in (0, 3, 4) for j in range(3) for k in range(2)}
    assert seeds(root) == expected, seeds(root)
    for r in root.findall('.//s:rect[@class="seed"][@data-target]', NS):
        bp = cells[int(r.attrib['data-target']), int(r.attrib['data-query'])]
        assert float(r.attrib['x']) == float(bp.attrib['x'])+1
        assert float(r.attrib['y']) == float(bp.attrib['y'])+1
    for r in frame:
        idx = int(r.attrib['data-index'])
        axis, bp = ('y', cells[idx, 1]) if r.attrib['data-strand'] == 'target' else ('x', cells[1, idx])
        assert float(r.attrib[axis]) == float(bp.attrib[axis])
    # Requesting SVG must not affect CSV or ensemble/ranked output.
    baseline = run(options=['--outNumber=1', '--outMode=C', '--outCsvCols=Zall,Eall'])
    plotted = run(options=['--outNumber=1', '--outMode=C', '--outCsvCols=Zall,Eall', f'--out=bpsvg:{svg}'])
    assert baseline.stdout == plotted.stdout
    # The green rectangle covers the reported MFE region, including when
    # ranked output is disabled. Query coordinates retain their original order.
    for target, query, options in (
        ('GGAGGG', 'CCCC', []),
        ('AAGGGGAA', 'AACCCCAA', []),
        ('GGAGGG', 'CCCC', ['--tIdxPos0=-4', '--qIdxPos0=-2']),
    ):
        reference = run(target, query, [*options, '--outNumber=1', '--outMode=C',
                                       '--outCsvCols=start1,end1,start2,end2'])
        row = next(csv.DictReader(reference.stdout.splitlines(), delimiter=';'))
        for count in (0, 1):
            mfe_plot = plot(target, query, [*options, f'--outNumber={count}'])
            outline = mfe_plot.find('.//s:rect[@class="mfe"][@data-target-start]', NS)
            assert outline is not None
            bounds = tuple(int(outline.attrib[k]) for k in (
                'data-target-start', 'data-target-end', 'data-query-start', 'data-query-end'))
            assert bounds == tuple(int(row[k]) for k in ('start1', 'end1', 'start2', 'end2'))
            ts, te, qs, qe = bounds
            matrix = pairs(mfe_plot)
            assert float(outline.attrib['x']) == float(matrix[ts, qs].attrib['x'])
            assert float(outline.attrib['y']) == float(matrix[te, qs].attrib['y'])
            assert float(outline.attrib['width']) == float(matrix[ts, qe].attrib['x'])-float(matrix[ts, qs].attrib['x'])+20
            assert float(outline.attrib['height']) == float(matrix[ts, qs].attrib['y'])-float(matrix[te, qs].attrib['y'])+20
            assert 'Green outline:' in ''.join(mfe_plot.itertext())
    unseeded = plot(opts=['--noSeed', f'--out=bpProb:{table}'])
    assert not seeds(unseeded)
    rows = list(csv.reader(table.open(), delimiter=';'))
    assert all(float(pairs(unseeded)[i,j].attrib['data-probability']) == float(rows[i][j])
               for i in range(1,7) for j in range(1,5))
    for seed_opts in ([], ['--noSeed'], ['--seedMaxUP=1']):
        approximate=plot(opts=['--mode=H', f'--out=bpProb:{table}', *seed_opts])
        assert 'Approximate' in approximate.find('s:title', NS).text
        rows=list(csv.reader(table.open(), delimiter=';'))
        assert rows[0][0]=='bpProbApproximate'
        assert all(float(pairs(approximate)[i,j].attrib['data-probability']) == float(rows[i][j])
                   for i in range(1,7) for j in range(1,5))
        if '--noSeed' in seed_opts:
            assert not seeds(approximate)
    # Explicit singleton/mixed seed masks exclude extension-only cells.
    opts = ['--seedTQ=1|&4|,2||&2||']
    mixed = plot('GGGG', 'CCCC', opts)
    assert seeds(mixed) == {(1,4), (2,3), (3,2)}
    restricted = plot('GGGG', 'CCCC', [*opts, '--tRegion=2-4', '--qRegion=1-3'])
    assert seeds(restricted) == {(2,3), (3,2)}
    # A seed crossing a searched region boundary is not annotated.
    cut = plot('GGGG', 'CCCC', ['--seedTQ=2||&2||', '--tRegion=1-2,3-4'])
    assert not seeds(cut)
    assert all(r.attrib['data-probability'] == 'NA' for r in pairs(cut).values())
    # Four disjoint rectangles merge both probabilities and masks.
    regional = plot('GGGG', 'CCCC', ['--tRegion=1-2,3-4', '--qRegion=1-2,3-4'])
    assert seeds(regional) == {(1,2),(2,1),(1,4),(2,3),(3,2),(4,1),(3,4),(4,3)}
    # Signed guides follow displayed coordinates; there is no zero nucleotide.
    signed = plot('A'*72, 'A'*18, ['--tIdxPos0=-15', '--qIdxPos0=-15'])
    signed_cells = pairs(signed)
    for strand, length in (('target',72), ('query',18)):
        guides = signed.findall(f'.//s:line[@data-strand="{strand}"]', NS)
        expected_indices = {-10, 1} | ({10,20,30,40,50} if length == 72 else set())
        assert {int(g.attrib['data-index']) for g in guides} == expected_indices
        assert not any(0 in k for k in signed_cells)
        for g in guides:
            idx = int(g.attrib['data-index'])
            r = signed_cells[idx,-15] if strand == 'target' else signed_cells[-15,idx]
            axis = 'y' if strand == 'target' else 'x'
            expected_pos = float(r.attrib[axis]) + (0 if idx > 1 else 20) if strand == 'target' else float(r.attrib[axis]) + (20 if idx > 1 else 0)
            assert float(g.attrib[axis+'1']) == expected_pos
            assert ('major' in g.attrib['class']) == (idx % 50 == 0)
            assert ('origin' in g.attrib['class']) == (idx == 1)
    negative = plot('A'*13, 'A'*13, ['--tIdxPos0=-55', '--qIdxPos0=45'])
    assert 'major' in negative.find('.//s:line[@data-strand="target"][@data-index="-50"]', NS).attrib['class']
    # User sequence identifiers must round-trip through XML escaping.
    name = '<target & "quoted">'
    named = plot(opts=[f'--tId={name}'])
    assert name in named.find('s:title', NS).text
    assert named.find('.//s:text[@class="subtitle"]', NS).text == name+' / query'
    assert all(name in r.find('s:title', NS).text for r in pairs(named).values())
    assert all(name in r.find('s:title', NS).text for r in named.findall(
        './/s:rect[@data-type="intramolecular-pairing"][@data-strand="target"]', NS))
    assert name+" (5' to 3')" in [e.text for e in named.findall('.//s:text', NS)]
    # Accessibility uses its own per-strand scale and original query order.
    pu_t, pu_q = d/'target.pu', d/'query.pu'

    def one_nt(path):
        return [float(line.split()[1]) for line in path.read_text().splitlines()
                if line.strip() and not line.lstrip().startswith('#')]

    def check_frame(root, target_pu, query_pu):
        for strand, values in (('target', target_pu), ('query', query_pu)):
            assert any(0 < v < .9 for v in values)
            for r in root.findall(f'.//s:rect[@data-type="intramolecular-pairing"][@data-strand="{strand}"]', NS):
                assert math.isclose(float(r.attrib['data-probability']), 1-values[int(r.attrib['data-index'])-1],
                                    rel_tol=2e-5, abs_tol=5e-6), (strand, r.attrib, values)

    for temperature in (22, 37):
        acc_options = [f'--temperature={temperature}', f'--out=tPu:{pu_t}', f'--out=qPu:{pu_q}']
        native = plot('GGGGAAAACCCC', 'GCGCGAAAACGCGC', ['--energy=V', '--acc=C', *acc_options])
        check_frame(native, one_nt(pu_t), one_nt(pu_q))
        input_t, input_q = d/'input-target.pu', d/'input-query.pu'
        input_t.write_bytes(pu_t.read_bytes()); input_q.write_bytes(pu_q.read_bytes())
        for t_mode, q_mode in (('P', 'P'), ('P', 'C'), ('C', 'P'), ('E', 'E'), ('C', 'C')):
            ed_t, ed_q = d/'target.ed', d/'query.ed'
            options = ['--energy=B', f'--tAcc={t_mode}', f'--qAcc={q_mode}', *acc_options]
            if t_mode in ('P', 'E'):
                options.append(f'--tAccFile={input_t if t_mode == "P" else ed_t}')
            if q_mode in ('P', 'E'):
                options.append(f'--qAccFile={input_q if q_mode == "P" else ed_q}')
            # The first imported-Pu run also supplies ED fixtures for the E case.
            if (t_mode, q_mode) in (('P', 'P'), ('C', 'C')):
                options += [f'--out=tAcc:{ed_t}', f'--out=qAcc:{ed_q}']
            converted = plot('GGGGAAAACCCC', 'GCGCGAAAACGCGC', options)
            check_frame(converted, one_nt(pu_t), one_nt(pu_q))
            for mode, input_pu, output_pu in ((t_mode, input_t, pu_t), (q_mode, input_q, pu_q)):
                if mode in ('P', 'E'):
                    # Independent input reference; ED quantization loses <2%
                    # at either temperature, but using interaction RT loses far more.
                    assert all(math.isclose(actual, original, rel_tol=.02, abs_tol=1e-10)
                               for actual, original in zip(one_nt(output_pu), one_nt(input_pu)))
            if t_mode == q_mode == 'C':
                # Nussinov accessibility has RT=1, independent of temperature.
                check_frame(converted, [math.exp(-ed) for ed in one_nt(ed_t)],
                            [math.exp(-ed) for ed in one_nt(ed_q)])
    for t,q in (('AAAA','AAAA'), ('G','C')):
        empty = plot(t,q)
        assert all(r.attrib['data-probability'] == 'NA' for r in pairs(empty).values())
        assert not seeds(empty)
        assert not empty.findall('.//s:rect[@class="mfe"]', NS)
        assert 'empty interaction ensemble' in ''.join(empty.itertext())
    compressed = d/'pairs.svg.gz'
    run(options=[f'--out=bpsvg:{compressed}'])
    assert len(pairs(load(compressed))) == 24
    # Validation/numerical failures must not publish a partial document.
    for opts in (['--mode=S'], ['--seedMaxUP=1'], ['--model=X', '--mode=K'],
                 ['--windowWidth=3'], ['--seedTQ=1|&4|','--energyAdd=900']):
        failure = d/'failure.svg'
        run(options=[f'--out=bpsvg:{failure}', *opts], success=False)
        assert not failure.exists()
    run(options=['--out=bpsvg:/dev/full'], success=False)
    run(options=[f'--out=bpsvg:{svg}', f'--out=bpProb:{svg}'], success=False)
    # Multi-FASTA creates one SVG per pair; shared streams emit indivisible roots.
    fasta=d/'targets.fa'; fasta.write_text('>a\nGGGG\n>b\nGAGG\n')
    for threads in (1,2):
        prefix=d/f'multi{threads}.svg'
        run(str(fasta), 'CCCC', [f'--threads={threads}', f'--out=bpsvg:{prefix}'])
        files=sorted(d.glob(f'multi{threads}*.svg'))
        assert len(files)==2
        assert all(len(pairs(load(file)))==16 for file in files)
        for stream in ('STDOUT', 'STDERR'):
            shared=run(str(fasta), 'CCCC', [f'--threads={threads}', f'--out=bpsvg:{stream}'])
            data=shared.stdout if stream=='STDOUT' else shared.stderr
            blocks=re.findall(r'<svg\b.*?</svg>', data, re.S)
            assert len(blocks)==2
            assert all(len(pairs(ET.fromstring(block)))==16 for block in blocks)

print('base-pair probability SVG regressions passed')

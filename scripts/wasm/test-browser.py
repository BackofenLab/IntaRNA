"""Exercise the staged browser interface with real WASM at a nested URL."""
import argparse
from pathlib import Path
from io import BytesIO
import re
import json
from urllib.parse import parse_qs, urljoin, urlparse
from zipfile import ZipFile
from playwright.sync_api import sync_playwright, expect

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('url')
parser.add_argument('--browser', choices=['chromium', 'firefox'], default='chromium')
parser.add_argument('--screenshot')
args = parser.parse_args()

with sync_playwright() as playwright:
    browser = getattr(playwright, args.browser).launch(headless=True)
    page = browser.new_page(viewport={'width': 1440, 'height': 1100})
    errors = []
    csp_errors = []
    page.on('pageerror', lambda error: errors.append(str(error)))
    page.on('console', lambda message: csp_errors.append(message.text) if message.type == 'error' and ('Content Security Policy' in message.text or 'Refused to apply inline style' in message.text) else None)
    # Keep browser tests offline and inspect the exact third-party URL separately.
    visualization_urls = []
    def viewer(route):
        visualization_urls.append(route.request.url)
        route.fulfill(content_type='text/html', body='<h1>vaRRI viewer</h1>')
    page.route('https://backofenlab.github.io/vaRRI/**', viewer)
    page.goto(args.url)
    page.wait_for_load_state('networkidle')
    compute = page.locator('#compute')
    status = page.locator('#result-status')
    sessions = page.locator('.session')
    active = page.locator('.session:has([data-role=panel]:not([hidden]))')
    rows = active.locator('tbody tr')

    def open_group(key):
        page.locator('#tab-' + key).click()
        expect(page.locator('#' + key)).to_be_visible()
        expect(page.locator('#tab-' + key)).to_have_attribute('aria-selected', 'true')
        assert page.locator('.parameter-group:visible').count() == 1

    def load_example(task='prediction'):
        page.locator('#examples-menu summary').click()
        page.get_by_role('button', name='Evaluate a duplex' if task == 'evaluation' else 'Short duplex', exact=False).click()
        expect(page.locator('#query')).to_have_value('>query\nGGGAGGAGGG\n')

    def download_archive(session=None):
        with page.expect_download() as downloaded:
            (active if session is None else session).locator('[data-role=download-zip]').click()
        with ZipFile(BytesIO(Path(downloaded.value.path()).read_bytes())) as archive:
            assert archive.testzip() is None
            return {name: archive.read(name).decode() for name in archive.namelist()}

    def fill(key, value):
        parent = page.locator('#' + key).evaluate('node => node.closest(".parameter-group")?.id')
        if parent:
            open_group(parent)
        page.locator('#' + key).fill(str(value))

    def predict():
        count = sessions.count()
        compute.click()
        expect(status).to_contain_text('complete:', timeout=120000)
        expect(page.locator('#input-error')).not_to_be_visible()
        expect(sessions).to_have_count(count + 1)

    expect(page.locator('#query')).to_have_value('')
    expect(page.locator('#target')).to_have_value('')
    expect(page.locator('#query')).to_have_attribute('placeholder', '>query\nGGGAGGAGGG')
    compute.click()
    expect(page.locator('#input-error')).to_contain_text('enter a sequence')
    metadata = page.request.get(urljoin(args.url, 'build-info.json')).json()
    vienna_version = metadata['viennarna']
    assert re.fullmatch(r'2\.\d+\.\d+', vienna_version)
    assert metadata['viennarnaArchive']['version'] == vienna_version
    expect(page.locator('#version-info')).to_contain_text(f'ViennaRNA {vienna_version}')
    for logo in page.locator('.site-header img').all():
        assert logo.evaluate('node => node.complete && node.naturalWidth > 0')
    for key in ['seed-options', 'output-options', 'csv-options', 'interaction-options', 'query-options', 'target-options', 'general-options', 'sequence-input']:
        open_group(key)

    expect(page.get_by_role('radio', name='Prediction')).to_be_checked()
    expect(page.get_by_role('radio', name='Evaluation')).to_be_visible()
    page.get_by_role('radio', name='Prediction').focus()
    page.keyboard.press('ArrowRight')
    expect(page.get_by_role('radio', name='Evaluation')).to_be_checked()
    expect(page.locator('#task-help')).to_contain_text('Evaluate supplied structures')
    page.keyboard.press('ArrowLeft')
    expect(page.get_by_role('radio', name='Prediction')).to_be_checked()
    page.locator('#examples-menu summary').click()
    expect(page.get_by_role('button', name='Short duplex', exact=False)).to_be_visible()
    expect(page.get_by_role('button', name='Evaluate a duplex', exact=False)).not_to_be_visible()
    page.get_by_role('radio', name='Evaluation').check()
    expect(page.get_by_role('button', name='Evaluate a duplex', exact=False)).to_be_visible()
    expect(page.get_by_role('button', name='Short duplex', exact=False)).not_to_be_visible()
    page.get_by_role('radio', name='Prediction').check()
    page.locator('#examples-menu summary').click()

    # Tabs support keyboard navigation and wrap without page overflow.
    page.locator('#tab-sequence-input').focus()
    page.keyboard.press('ArrowRight')
    expect(page.locator('#tab-seed-options')).to_be_focused()
    expect(page.locator('#seed-options')).to_be_visible()
    page.keyboard.press('End')
    expect(page.locator('#tab-general-options')).to_be_focused()
    page.keyboard.press('Home')
    expect(page.locator('#tab-sequence-input')).to_be_focused()

    # Review layout, task defaults and modes work with native parameters.
    assert page.locator('[role=tab]').all_text_contents() == ['Sequences', 'Seed', 'Output', 'CSV', 'Model & Energy', 'Query', 'Target', 'CLI-only']
    expect(page.locator('.brand')).to_contain_text('IntaRNA package')
    expect(page.locator('#results-title')).to_have_text('Results')
    expect(page.locator('#outNumber')).to_have_value('10')
    expect(page.locator('#outNumberPerSeed')).to_have_value('0')
    assert page.locator('.compute-bar').evaluate('node => node.firstElementChild.id') == 'compute'
    assert page.locator('.sequence-meta').first.evaluate('node => node.firstElementChild.tagName') == 'LABEL'
    fill('query', '>named\tignored description\n' + 'G' * 1000 + '\n' * 20 + 'GG')
    expect(page.locator('#query-warning')).to_contain_text('ignoring FASTA description')
    expect(page.locator('#query-warning')).to_contain_text('named')
    assert page.locator('#query').evaluate('node => node.scrollWidth > node.clientWidth && node.scrollHeight > node.clientHeight')
    fill('query', '>named   \nGGGAGGAGGG')
    expect(page.locator('#query-warning')).not_to_be_visible()

    for prefix, group in [('q', 'query-options'), ('t', 'target-options')]:
        open_group(group)
        # Each explicit row shares label height, including wrapped labels, and
        # compact controls stay left-aligned instead of filling the panel.
        for keys in [[prefix + suffix for suffix in ['Set', 'Region', 'IntLenMax', 'IntLoopMax', 'RegionLenMax']], [prefix + suffix for suffix in ['Acc', 'AccW', 'AccL', 'PfScale']]]:
            tops = [page.locator('#' + key).bounding_box()['y'] for key in keys]
            assert max(tops) - min(tops) < 1, tops
            label_bottoms = [page.locator('#' + key).locator('xpath=preceding-sibling::label').bounding_box()['y'] + page.locator('#' + key).locator('xpath=preceding-sibling::label').bounding_box()['height'] for key in keys]
            assert max(label_bottoms) - min(label_bottoms) < 1, label_bottoms
        first = page.locator('#' + prefix + 'Set').bounding_box()
        second = page.locator('#' + prefix + 'Region').bounding_box()
        assert 8 <= second['x'] - first['x'] - first['width'] <= 20
        assert page.locator('#' + prefix + 'AccW').bounding_box()['width'] <= 224
        expect(page.locator('#' + prefix + 'AccFile')).to_be_disabled()
        page.locator('#' + prefix + 'Acc').select_option('P')
        expect(page.locator('#' + prefix + 'AccFile')).to_be_enabled()
        for suffix in ['AccW', 'AccL']:
            expect(page.locator('#' + prefix + suffix)).to_be_disabled()
        page.locator('#' + prefix + 'Acc').select_option('N')
        for suffix in ['AccW', 'AccL', 'AccConstr', 'AccFile', 'PfScale', 'Shape', 'ShapeMethod', 'ShapeConversion']:
            expect(page.locator('#' + prefix + suffix)).to_be_disabled()
        expect(page.locator('input[data-for="' + prefix + 'Shape"]')).to_be_disabled()
        page.locator('#' + prefix + 'Acc').select_option('C')
        expect(page.locator('#' + prefix + 'Shape')).to_be_enabled()
    page.locator('#accNoLP-target').check()
    expect(page.locator('#accNoLP')).to_be_checked()
    open_group('query-options')
    page.locator('#accNoLP').uncheck()
    expect(page.locator('#accNoLP-target')).not_to_be_checked()
    open_group('interaction-options')
    expect(page.locator('#helixMinBP')).to_be_disabled()
    page.locator('#model').select_option('B')
    expect(page.locator('#helixMinBP')).to_be_enabled()
    page.locator('#model').select_option('X')
    expect(page.locator('#helixMinBP')).to_be_disabled()
    open_group('seed-options')
    expect(page.locator('#seedTQ')).to_be_disabled()
    seed_box = page.locator('#seed-mode').bounding_box()
    explicit_box = page.locator('#seedTQ').bounding_box()
    assert abs(seed_box['y'] - explicit_box['y']) < 5
    assert seed_box['x'] + seed_box['width'] < explicit_box['x']
    assert seed_box['width'] < explicit_box['width']
    page.locator('#seed-mode').select_option('none')
    for key in ['seedBP', 'seedTQ', 'seedQRange', 'seedNoGU']:
        expect(page.locator('#' + key)).to_be_disabled()
    page.locator('#seed-mode').select_option('explicit')
    expect(page.locator('#seedTQ')).to_be_enabled()
    expect(page.locator('#seedBP')).to_be_disabled()
    page.locator('#seed-mode').select_option('identify')
    expect(page.locator('#seedBP')).to_be_enabled()
    expect(page.locator('#seedTQ')).to_be_disabled()
    open_group('output-options')
    output_rows = [
        [page.locator('#' + key).locator('xpath=..').bounding_box() for key in keys]
        for keys in [
            ['outNumber', 'outNumberPerSeed', 'outOverlap'],
            ['outMaxE', 'outDeltaE', 'outMinPu'],
            ['outBestSeedOnly', 'outNoLP', 'outNoGUend', 'outPerRegion', 'outPairwise'],
        ]
    ]
    for boxes in output_rows:
        assert max(box['y'] for box in boxes) - min(box['y'] for box in boxes) < 1, boxes
    assert output_rows[0][0]['x'] < output_rows[0][1]['x'] < output_rows[0][2]['x']
    assert output_rows[0][0]['y'] < output_rows[1][0]['y'] < output_rows[2][0]['y']
    output_labels = page.locator('#output-options .checkbox-grid label')
    assert output_labels.count() == 15
    assert all(label.get_attribute('title') for label in output_labels.all())
    expect(page.locator('#output-options input[value="pMinE"]').locator('xpath=..')).to_have_attribute('title', re.compile('energy heatmap'))
    open_group('csv-options')
    assert page.locator('#outCsvSort').bounding_box()['width'] < 310
    groups = page.locator('.csv-group')
    assert groups.count() == 5
    for group in groups.all():
        group.get_by_role('button', name=re.compile('^Select all')).click()
        assert group.locator('input:not(:checked)').count() == 0
        group.get_by_role('button', name='Deselect all').click()
        assert group.locator('input:checked').count() == 0
    page.get_by_role('button', name='Task defaults', exact=True).click()
    for key in ['ED1', 'ED2', 'seedE', 'seedStart1', 'seedStart2']:
        expect(page.locator('.csv-group input[value="' + key + '"]')).to_be_checked()
    assert all(label.get_attribute('title') and ' — ' not in label.inner_text() for label in groups.locator('label').all())
    # Inspect all changed panels at desktop and mobile widths.
    for width in [1440, 390]:
        page.set_viewport_size({'width': width, 'height': 1100 if width == 1440 else 844})
        for key in ['sequence-input', 'seed-options', 'output-options', 'csv-options', 'interaction-options', 'query-options', 'target-options', 'general-options']:
            open_group(key)
            assert page.evaluate('document.documentElement.scrollWidth <= innerWidth'), (width, key)
            if args.screenshot:
                path = Path(args.screenshot)
                page.screenshot(path=str(path.with_name(f'{path.stem}-{key}-{width}.png')), full_page=True)
    page.set_viewport_size({'width': 1440, 'height': 1100})
    page.locator('#reset').click()

    # Upload and file drop use local text, validate before replacing, and reset.
    with page.expect_file_chooser() as chooser:
        page.locator('#query-upload-button').click()
    chooser.value.set_files({'name': 'query.fa', 'mimeType': 'text/plain', 'buffer': b'>uploaded\nGGGAGGAGGG\n'})
    expect(page.locator('#query')).to_have_value('>uploaded\nGGGAGGAGGG\n')
    page.locator('#target').evaluate(r"""node => {
      const transfer = new DataTransfer();
      transfer.items.add(new File(['>dropped\nCCCUCCUCCC\n'], 'target.fa', {type:'text/plain'}));
      node.dispatchEvent(new DragEvent('drop', {dataTransfer:transfer, bubbles:true, cancelable:true}));
    }""")
    expect(page.locator('#target')).to_have_value('>dropped\nCCCUCCUCCC\n')
    page.locator('#query-upload').set_input_files({'name': 'bad.fa', 'mimeType': 'text/plain', 'buffer': b'not a sequence!'})
    expect(page.locator('#input-error')).to_be_visible()
    expect(page.locator('#query')).to_have_value('>uploaded\nGGGAGGAGGG\n')
    # A file read completing after Reset cannot restore the previous input.
    page.evaluate(r"""() => { window.originalFileText = File.prototype.text; File.prototype.text = function() { return new Promise(resolve => { window.finishRead = () => resolve('>late\nGGGAGGAGGG\n'); }); }; }""")
    page.locator('#query-upload').set_input_files({'name': 'late.fa', 'mimeType': 'text/plain', 'buffer': b'>late\nGGGAGGAGGG\n'})
    page.locator('#reset').click()
    page.evaluate('() => { File.prototype.text = window.originalFileText; window.finishRead(); }')
    expect(page.locator('#query')).to_have_value('')

    # JavaScript validation works even when the offending control is in another tab.
    fill('query', '>bad\nACGU-ACGU')
    compute.click()
    expect(page.locator('#input-error')).to_contain_text('use only')
    expect(page.locator('#loading')).not_to_be_visible()
    load_example()
    fill('temperature', '101')
    open_group('sequence-input')
    compute.click()
    expect(page.locator('#input-error')).to_contain_text('Temperature')
    page.locator('#reset').click()
    expect(page.locator('#query')).to_have_value('')
    expect(page.locator('#target')).to_have_value('')
    expect(page.locator('#rri')).to_have_value('')
    load_example()

    # Disabled stale edits neither block a run nor leak into JSON. Re-enabling
    # restores local edits; loading/storing uses effective native settings.
    fill('seedBP', '')
    page.locator('#seed-mode').select_option('explicit')
    compute.click()
    expect(page.locator('#input-error')).to_contain_text('Explicit seed: enter')
    fill('seedTQ', '1|||&1|||')
    predict()
    explicit_settings = json.loads(download_archive()['parameter.json'])
    assert explicit_settings['parameters']['seedTQ'] == '1|||&1|||'
    assert 'seedBP' not in explicit_settings['parameters']
    with page.expect_download() as stored:
        page.locator('#store-parameters').click()
    explicit_file = Path(stored.value.path()).read_bytes()
    page.locator('#reset').click()
    page.locator('#parameter-upload').set_input_files({'name': 'explicit.json', 'mimeType': 'application/json', 'buffer': explicit_file})
    expect(page.locator('#seed-mode')).to_have_value('explicit')
    open_group('seed-options')
    page.locator('#seedTQ').fill('invalid inactive seed')
    page.locator('#seed-mode').select_option('none')
    fill('qShape', 'invalid inactive shape')
    page.locator('#qAcc').select_option('N')
    open_group('target-options')
    page.locator('#accNoGUend-target').check()
    predict()
    inactive_settings = json.loads(download_archive()['parameter.json'])['parameters']
    assert inactive_settings['noSeed'] is True
    assert inactive_settings['qAcc'] == 'N'
    assert inactive_settings['accNoGUend'] is True
    assert 'qShape' not in inactive_settings and 'seedTQ' not in inactive_settings
    open_group('query-options')
    page.locator('#qAcc').select_option('C')
    expect(page.locator('#qShape')).to_have_value('invalid inactive shape')
    page.locator('#reset').click()
    load_example()

    # Reference prediction, selectable CSV output and independent table visibility.
    fill('outNumber', 10)
    open_group('csv-options')
    page.locator('#outCsvSort').select_option('E')
    predict()
    expect(rows).to_have_count(10)
    assert active.locator('.table-scroll').evaluate('''node => {
        const body = node.querySelector('tbody');
        const copies = [...body.children].map(row => row.cloneNode(true));
        body.append(...copies);
        const scrolls = node.clientHeight <= 500 && node.scrollHeight > node.clientHeight
            && getComputedStyle(node).overflowY === 'auto';
        copies.forEach(copy => copy.remove());
        return scrolls;
    }''')
    assert rows.first.locator('td').all_text_contents()[1:7] == ['target', '1', '10', 'query', '1', '10']
    assert rows.first.locator('td').nth(9).inner_text() == '-21.40'
    assert rows.first.locator('td').last.inner_text() == '0.00'
    page.get_by_role('button', name='Sort by Energy (kcal/mol)', exact=False).click()
    assert float(rows.first.locator('td').nth(9).inner_text()) == -14.3
    page.get_by_role('button', name='Sort by Energy (kcal/mol)', exact=False).click()
    active.locator('[data-role=search]').fill('-21.4')
    expect(rows).to_have_count(1)
    active.locator('[data-role=search]').fill('no match')
    expect(active.locator('[data-role=no-matches]')).to_be_visible()
    active.locator('[data-role=search]').fill('')
    active.locator('.column-picker summary').click()
    active.locator('[data-role=visible-columns] input[value="start1"]').uncheck()
    assert 'Target start' not in active.locator('thead').inner_text()
    active.locator('[data-role=visible-columns] input[value="start1"]').check()
    with page.expect_download() as downloaded:
        active.locator('[data-role=download]').click()
    csv = Path(downloaded.value.path()).read_text()
    assert csv.splitlines()[0] == 'id1;start1;end1;id2;start2;end2;subseqDP;hybridDP;E;ED1;ED2;seedE;seedStart1;seedStart2'
    assert 'target;1;10;query;1;10;' in csv
    assert 'deltaE' not in csv
    assert active.locator('thead th').first.inner_text() == ''
    assert active.locator('thead th button br').count() >= 4
    assert rows.first.locator('td').first.evaluate('node => getComputedStyle(node).paddingLeft') == '5px'
    assert not visualization_urls  # Third-party content loads only on selection.
    rows.first.get_by_role('button', name='Visualize').click()
    expect(active.locator('[data-role=varri-frame]')).to_be_visible()
    expect(active.frame_locator('[data-role=varri-frame]').get_by_role('heading', name='vaRRI viewer')).to_be_visible()
    expect(rows.first.locator('button')).to_have_attribute('aria-pressed', 'true')
    query = parse_qs(urlparse(active.locator('[data-role=varri-frame]').get_attribute('src')).query)
    assert query['sequence'] == ['CCCUCCUCCC&GGGAGGAGGG']
    assert query['structure'] == ['((((((((((&))))))))))']
    assert query['seqName1'] == ['target'] and query['seqName2'] == ['query']
    assert query['startIndex1'] == ['1']
    assert query['showRenderingOnly'] == ['true']
    assert 'showRenderingOnly' not in active.locator('[data-role=varri-link]').get_attribute('href')
    archive = download_archive()
    assert archive['intarna-interactions.csv'] == csv
    saved = json.loads(archive['parameter.json'])
    assert set(archive) == {'parameter.json', 'intarna-interactions.csv'}
    assert saved['parameters']['query']['text'] == '>query\nGGGAGGAGGG\n'
    assert 'IntaRNA 3.4.1' in archive['parameter.json']
    assert f'Vienna RNA package {vienna_version}' in archive['parameter.json']
    assert saved['parameters']['outNumber'] == 10
    assert 'IntaRNA --parameterFile=parameter.json' in archive['parameter.json']
    assert 'temperature' not in saved['parameters']
    # Store snapshots the editable form; Load is atomic and preserves CSV order.
    with page.expect_download() as stored:
        page.locator('#store-parameters').click()
    settings = Path(stored.value.path()).read_bytes()
    imported = json.loads(settings)
    imported['parameters']['outCsvCols'] = 'E,id2'
    imported['parameters']['temperature'] = 24
    page.locator('#parameter-upload').set_input_files({'name': 'settings.json', 'mimeType': 'application/json', 'buffer': json.dumps(imported).encode()})
    expect(page.locator('#temperature')).to_have_value('24')
    expect(active.locator('[data-role=results-view]')).to_be_visible()
    with page.expect_download() as stored_again:
        page.locator('#store-parameters').click()
    assert json.loads(Path(stored_again.value.path()).read_text())['parameters']['outCsvCols'] == 'E,id2'
    page.locator('#parameter-upload').set_input_files({'name': 'bad.json', 'mimeType': 'application/json', 'buffer': b'{invalid'})
    expect(page.locator('#input-error')).to_be_visible()
    expect(page.locator('#temperature')).to_have_value('24')
    page.locator('#parameter-upload').set_input_files({'name': 'settings.json', 'mimeType': 'application/json', 'buffer': settings})
    expect(page.locator('#temperature')).to_have_value('37')
    predict()
    if args.screenshot:

        open_group('sequence-input')
        page.screenshot(path=args.screenshot, full_page=True)
    page.set_viewport_size({'width': 390, 'height': 844})
    assert page.evaluate('document.documentElement.scrollWidth <= innerWidth')
    page.set_viewport_size({'width': 1440, 'height': 1100})

    # Multi-FASTA prediction, signed inclusive output and additional files.
    load_example()
    fill('outNumber', 1)
    fill('query', '>q1\nGGGAGGAGGG\n>q2\nGGGAGGAGGG')
    fill('qIdxPos0', -3)
    fill('tIdxPos0', -5)
    open_group('output-options')
    page.locator('#output-options input[value="qMinE"]').check()
    predict()
    expect(rows).to_have_count(2)
    assert rows.first.locator('td').all_text_contents()[1:7] == ['target', '-5', '5', 'q1', '-3', '7']
    output_files = active.locator('[data-role=file-list] details')
    expect(output_files).to_have_count(2)
    expect(output_files.nth(0).locator('summary')).to_contain_text('target : q1')
    expect(output_files.nth(1).locator('summary')).to_contain_text('target : q2')
    output_files.first.locator('summary').click()
    expect(output_files.first.get_by_role('tab', name='Visualization')).to_have_attribute('aria-selected', 'true')
    assert output_files.first.get_by_role('tab', name='Visualization').evaluate('node => getComputedStyle(node).backgroundColor') == page.locator('#tab-output-options').evaluate('node => getComputedStyle(node).backgroundColor')
    expect(output_files.first.locator('.js-plotly-plot')).to_be_visible(timeout=30000)
    expect(output_files.first.locator('.gtitle')).to_contain_text('q1')
    if args.screenshot:
        path = Path(args.screenshot)
        page.screenshot(path=str(path.with_name(f'{path.stem}-additional-output.png')), full_page=True)
    with page.expect_download(timeout=30000) as image:
        output_files.first.get_by_role('button', name='Export PNG').click()
    assert Path(image.value.path()).read_bytes().startswith(b'\x89PNG\r\n\x1a\n')
    output_files.first.get_by_role('tab', name='File Content').click()
    expect(output_files.first.locator('pre')).to_be_visible()
    expect(output_files.first.locator('pre')).to_contain_text('minE')
    with page.expect_download() as output_download:
        output_files.first.get_by_role('button', name=re.compile('^Download qMinE')).click()
    assert 'minE' in Path(output_download.value.path()).read_text()
    archive = download_archive()
    assert len([name for name in archive if name.startswith('qMinE')]) == 2
    saved = json.loads(archive['parameter.json'])
    assert saved['parameters']['query']['text'].count('>') == 2
    assert saved['parameters']['qIdxPos0'] == -3
    rows.nth(1).get_by_role('button', name='Visualize').click()
    query = parse_qs(urlparse(active.locator('[data-role=varri-frame]').get_attribute('src')).query)
    assert query['seqName2'] == ['q2']
    assert query['startIndex1'] == ['-5'] and query['startIndex2'] == ['-3']

    # FASTA identity follows native first-token naming, including collisions.
    fill('query', '>q first\nGGGAGGAGGG\n>q second\nGGGAGGAGGG')
    compute.click()
    expect(page.locator('#input-error')).to_contain_text('unique')
    expect(active.locator('[data-role=varri-frame]')).to_be_visible()
    expect(active.locator('[data-role=download-zip]')).to_be_visible()
    assert download_archive() == archive

    # Evaluation resets native defaults and removes prediction-only controls.
    fill('query', '>q1\nGGGAGGAGGG\n>q2\nGGGAGGAGGG')
    page.get_by_role('radio', name='Evaluation').check()
    expect(page.locator('#seed-options')).not_to_be_visible()
    expect(page.locator('#outNumber')).to_be_disabled()
    fill('rri', '1||||||||||&1||||||||||')
    compute.click()
    expect(page.locator('#input-error')).to_contain_text('exactly one')
    load_example('evaluation')
    fill('rri', '1||||||||||&1||||||||||:1|||||||||&2|||||||||')
    open_group('output-options')
    page.locator('#output-options input[value="detailed"]').check()
    open_group('csv-options')
    page.locator('#csv-options input[value="ED1"]').check()
    predict()
    expect(status).to_contain_text('Evaluation complete')
    expect(rows).to_have_count(2)
    expect(active.locator('thead')).to_contain_text('ED1')
    expect(active.locator('[data-role=file-list] summary')).to_contain_text('intarna-detailed.txt')
    active.locator('[data-role=file-list] summary').click()
    expect(active.locator('[data-role=file-list] pre')).to_contain_text('-21.4')
    with page.expect_download() as report:
        active.locator('[data-role=file-list] button').click()
    assert '-21.4' in Path(report.value.path()).read_text()
    archive = download_archive()
    assert '-21.4' in archive['intarna-detailed.txt']
    saved = json.loads(archive['parameter.json'])
    assert saved['parameters']['personality'] == 'IntaRNAeval'
    assert saved['parameters']['rri'] == '1||||||||||&1||||||||||:1|||||||||&2|||||||||'
    assert 'seedBP' not in saved['parameters']
    assert 'detailed' in saved['gui']['reports']

    # Visualization stays available when structure/identity columns are deselected.
    open_group('csv-options')
    page.get_by_role('button', name='Clear selection', exact=True).click()
    page.locator('#csv-options input[value="E"]').check()
    predict()
    rows.first.get_by_role('button', name='Visualize').click()
    expect(active.locator('[data-role=varri-frame]')).to_be_visible()
    archive = download_archive()
    assert archive['intarna-interactions.csv'].splitlines()[0] == 'E'
    assert json.loads(archive['parameter.json'])['parameters']['outCsvCols'] == 'E'
    fill('rri', '1|||&1||')
    compute.click()
    expect(page.locator('#input-error')).to_contain_text('same nonzero')

    # Startup and synchronous C++ cancellation both allow a fresh prediction.
    page.get_by_role('radio', name='Prediction').check()
    load_example()
    held_routes = []
    session_count = sessions.count()
    page.route('**/wasm/intarna.wasm', lambda route: held_routes.append(route))
    compute.click()
    expect(page.locator('#loading')).to_be_visible()
    page.locator('#abort').click()
    expect(status).to_contain_text('aborted')
    expect(sessions).to_have_count(session_count)
    for route in held_routes:
        route.abort()
    page.unroute('**/wasm/intarna.wasm')
    fill('query', '>long_query\n' + 'G' * 1000)
    fill('target', '>long_target\n' + 'C' * 1000)
    compute.click()
    expect(page.locator('#loading-status')).to_contain_text('Predicting RNA', timeout=120000)
    page.keyboard.press('Escape')
    expect(status).to_contain_text('aborted')
    expect(sessions).to_have_count(session_count)
    load_example()
    predict()

    fill('query', '>no_pairs\nNNNNNNNNNN')
    compute.click()
    expect(status).to_contain_text('No favorable interactions', timeout=120000)
    expect(rows).to_have_count(0)
    expect(active.locator('[data-role=varri-frame]')).not_to_be_visible()
    assert download_archive()['intarna-interactions.csv'].count('\n') == 1
    fill('query', '><img/src=x/onerror=alert(1)>\nGGGAGGAGGG')
    predict()
    assert active.locator('tbody img').count() == 0
    assert '<img' in rows.first.inner_text()

    load_example()
    session_count = sessions.count()
    page.route('**/wasm/intarna.js', lambda route: route.abort())
    compute.click()
    expect(status).to_contain_text('failed', timeout=30000)
    expect(sessions).to_have_count(session_count)
    expect(page.locator('#loading')).not_to_be_visible()
    expect(compute).to_be_enabled()
    page.unroute('**/wasm/intarna.js')
    predict()
    # Sessions retain independent JSON inputs, outputs and exploration state.
    count_before_reload = sessions.count()
    assert count_before_reload > 1
    page.reload()
    page.wait_for_load_state('networkidle')
    expect(sessions).to_have_count(0)
    load_example()
    fill('temperature', 24)
    fill('qIdxPos0', -3)
    fill('tIdxPos0', -5)
    fill('qShape', '1 0.2\n2 0.3')
    fill('outNumber', 3)
    open_group('output-options')
    page.locator('#output-options input[value="qMinE"]').check()
    predict()
    first_id = active.get_attribute('data-session-id')
    first = page.locator(f'[data-session-id="{first_id}"]')
    first.get_by_label('Session name', exact=True).fill('Cool duplex')
    first.get_by_label('Session name', exact=True).press('Enter')
    first_archive = download_archive(first)
    first_csv = first_archive['intarna-interactions.csv']
    first.locator('[data-role=search]').fill('no match')
    expect(first.locator('[data-role=no-matches]')).to_be_visible()
    first.locator('.column-picker summary').click()
    first.locator('[data-role=visible-columns] input[value="start1"]').uncheck()
    first_diagnostics = first.locator('[data-role=diagnostic-text]').text_content()

    # Edits and Reset preserve outputs. Restore also works for an unfolded run.
    fill('temperature', 31)
    expect(first.locator('[data-role=panel]')).to_be_visible()
    assert download_archive(first) == first_archive
    first.locator('[data-role=summary]').click()
    expect(page.locator('#temperature')).to_have_value('24')
    fill('temperature', 31)
    first.locator('[data-role=toggle]').click()
    expect(page.locator('#temperature')).to_have_value('24')
    expect(first.locator('[data-role=panel]')).to_be_visible()
    first.locator('[data-role=toggle]').click()
    expect(first.locator('[data-role=panel]')).not_to_be_visible()
    first.locator('[data-role=toggle]').press('Enter')
    expect(first.locator('[data-role=panel]')).to_be_visible()
    page.locator('#reset').click()
    expect(sessions).to_have_count(1)
    expect(page.locator('#query')).to_have_value('')
    assert download_archive(first) == first_archive

    # A different task, CSV order and output report own a second result area.
    page.get_by_role('radio', name='Evaluation').check()
    load_example('evaluation')
    with page.expect_download() as stored:
        page.locator('#store-parameters').click()
    second_settings = json.loads(Path(stored.value.path()).read_text())
    second_settings['parameters']['temperature'] = 42
    second_settings['parameters']['outCsvCols'] = 'E,id2'
    second_settings['gui'] = {'reports': ['detailed']}
    page.locator('#parameter-upload').set_input_files({'name': 'evaluation.json', 'mimeType': 'application/json', 'buffer': json.dumps(second_settings).encode()})
    predict()
    expect(sessions).to_have_count(2)
    second_id = active.get_attribute('data-session-id')
    assert second_id != first_id
    second = page.locator(f'[data-session-id="{second_id}"]')
    expect(first.locator('[data-role=panel]')).not_to_be_visible()
    expect(first.locator('[data-role=download-zip]')).to_be_visible()
    assert download_archive(first) == first_archive
    expect(page.get_by_role('radio', name='Evaluation')).to_be_checked()
    expect(page.locator('#temperature')).to_have_value('42')
    second_archive = download_archive(second)
    assert second_archive['intarna-interactions.csv'].splitlines()[0] == 'E;id2'
    assert 'intarna-detailed.txt' in second_archive

    # Selecting reloads all computed inputs while retaining each table's filter.
    first.locator('.session-header').focus()
    first.locator('.session-header').press('Space')
    expect(second.locator('[data-role=panel]')).not_to_be_visible()
    expect(first.locator('[data-role=panel]')).to_be_visible()
    assert first.locator('.session-header').evaluate('node => node.firstElementChild.dataset.role') == 'toggle'
    expect(first.locator('[data-role=toggle]')).to_have_text('Hide')
    expect(page.get_by_role('radio', name='Prediction')).to_be_checked()
    expect(page.locator('#qIdxPos0')).to_have_value('-3')
    expect(page.locator('#tIdxPos0')).to_have_value('-5')
    expect(page.locator('#qShape')).to_have_value('1 0.2\n2 0.3')
    expect(page.locator('#temperature')).to_have_value('24')
    expect(page.locator('#output-options input[value="qMinE"]')).to_be_checked()
    expect(page.locator('#output-options input[value="detailed"]')).not_to_be_checked()
    expect(first.locator('[data-role=search]')).to_have_value('no match')
    expect(first.locator('[data-role=visible-columns] input[value="start1"]')).not_to_be_checked()
    assert first.locator('[data-role=diagnostic-text]').text_content() == first_diagnostics
    with page.expect_download() as restored:
        page.locator('#store-parameters').click()
    assert Path(restored.value.path()).read_text() == first_archive['parameter.json']
    first.locator('[data-role=search]').fill('')
    expect(first.locator('tbody tr')).to_have_count(3)
    first.locator('tbody tr').first.get_by_role('button', name='Visualize').click()
    first_url = first.locator('[data-role=varri-frame]').get_attribute('src')
    assert parse_qs(urlparse(first_url).query)['startIndex1'] == ['-5']

    # A pending file read must not overwrite a selected session's snapshot.
    page.evaluate(r"""() => { window.originalFileText = File.prototype.text; File.prototype.text = function() { return new Promise(resolve => { window.finishRead = () => resolve('>late\nGGGAGGAGGG\n'); }); }; }""")
    page.locator('#query-upload').set_input_files({'name': 'late.fa', 'mimeType': 'text/plain', 'buffer': b'>late\nGGGAGGAGGG\n'})
    second.locator('[data-role=toggle]').click()
    page.evaluate('() => { File.prototype.text = window.originalFileText; window.finishRead(); }')
    expect(page.locator('#query')).to_have_value('>query\nGGGAGGAGGG\n')
    expect(page.locator('#qShape')).to_have_value('')
    expect(page.locator('#temperature')).to_have_value('42')
    expect(page.get_by_role('radio', name='Evaluation')).to_be_checked()
    expect(page.locator('#rri')).to_have_value(second_settings['parameters']['rri'])
    expect(page.locator('#output-options input[value="detailed"]')).to_be_checked()
    assert first.locator('[data-role=varri-frame]').get_attribute('src') is None
    with page.expect_download() as restored:
        page.locator('#store-parameters').click()
    assert Path(restored.value.path()).read_text() == second_archive['parameter.json']
    assert download_archive(first)['intarna-interactions.csv'] == first_csv

    # Names are literal text and cannot introduce markup, duplicate IDs or paths.
    first.get_by_label('Session name', exact=True).fill('<img src=x onerror=alert(1)> / β')
    first.get_by_label('Session name', exact=True).press('Enter')
    expect(page.locator('#temperature')).to_have_value('24')
    expect(first.locator('[data-role=panel]')).to_be_visible()
    assert first.locator('img').count() == 0
    with page.expect_download() as renamed:
        first.locator('[data-role=download-zip]').click()
    assert '/' not in renamed.value.suggested_filename
    assert '\\' not in renamed.value.suggested_filename
    ids = page.locator('[id]').evaluate_all('nodes => nodes.map(node => node.id)')
    assert len(ids) == len(set(ids))
    first.get_by_label('Session name', exact=True).fill('')
    first.get_by_label('Session name', exact=True).press('Enter')
    expect(first.get_by_label('Session name', exact=True)).to_have_value('Prediction 1')
    first.get_by_label('Session name', exact=True).fill('Cool duplex — 24 °C')
    first.get_by_label('Session name', exact=True).press('Enter')
    if args.screenshot:
        page.screenshot(path=args.screenshot, full_page=True)
    page.set_viewport_size({'width': 390, 'height': 844})
    assert page.evaluate('document.documentElement.scrollWidth <= innerWidth')
    for session in [first, second]:
        expect(session.locator('[data-role=download-zip]')).to_be_visible()
        expect(session.locator('[data-role=delete]')).to_be_visible()
    page.set_viewport_size({'width': 1440, 'height': 1100})

    # Delete releases only that output; current inputs and other sessions survive.
    second.locator('[data-role=toggle]').click()
    first.locator('[data-role=delete]').click()
    expect(sessions).to_have_count(1)
    expect(page.locator('#temperature')).to_have_value('42')
    assert download_archive(second) == second_archive
    second.locator('[data-role=delete]').click()
    expect(sessions).to_have_count(0)
    expect(compute).to_be_focused()
    expect(page.locator('#temperature')).to_have_value('42')
    predict()
    assert active.get_attribute('data-session-id') not in [first_id, second_id]
    assert download_archive() == second_archive
    page.reload()
    page.wait_for_load_state('networkidle')
    expect(sessions).to_have_count(0)
    assert errors == [], errors
    assert csp_errors == [], csp_errors
    browser.close()
print(f'{args.browser}: real WASM, groups/defaults/versions, prediction/evaluation, FASTA/coordinates, tabs/keyboard, visualization URLs/selection, ZIP replay contents, columns/reports, mobile, abort/recovery and retained session lifecycle passed')

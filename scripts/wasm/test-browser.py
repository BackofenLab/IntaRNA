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
    page.on('pageerror', lambda error: errors.append(str(error)))
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
    rows = page.locator('tbody tr')

    def open_group(key):
        page.locator('#tab-' + key).click()
        expect(page.locator('#' + key)).to_be_visible()
        expect(page.locator('#tab-' + key)).to_have_attribute('aria-selected', 'true')
        assert page.locator('.parameter-group:visible').count() == 1

    def load_example(task='prediction'):
        page.locator('#examples-menu summary').click()
        page.get_by_role('button', name='Evaluate a duplex' if task == 'evaluation' else 'Short duplex', exact=False).click()
        expect(page.locator('#query')).to_have_value('>query\nGGGAGGAGGG\n')

    def download_archive():
        with page.expect_download() as downloaded:
            page.locator('#download-zip').click()
        with ZipFile(BytesIO(Path(downloaded.value.path()).read_bytes())) as archive:
            assert archive.testzip() is None
            return {name: archive.read(name).decode() for name in archive.namelist()}

    def fill(key, value):
        parent = page.locator('#' + key).evaluate('node => node.closest(".parameter-group")?.id')
        if parent:
            open_group(parent)
        page.locator('#' + key).fill(str(value))

    def predict():
        compute.click()
        expect(status).to_contain_text('complete:', timeout=120000)
        expect(page.locator('#input-error')).not_to_be_visible()

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
    for key in ['query-options', 'target-options', 'seed-options', 'helix-options', 'shape-options', 'interaction-options', 'output-options', 'general-options', 'sequence-input']:
        open_group(key)

    expect(page.get_by_role('radio', name='Prediction')).to_be_checked()
    expect(page.get_by_role('radio', name='Evaluation')).to_be_visible()
    page.get_by_role('radio', name='Prediction').focus()
    page.keyboard.press('ArrowRight')
    expect(page.get_by_role('radio', name='Evaluation')).to_be_checked()
    expect(page.locator('#task-help')).to_contain_text('Evaluate supplied structures')
    page.keyboard.press('ArrowLeft')
    expect(page.get_by_role('radio', name='Prediction')).to_be_checked()

    # Tabs support keyboard navigation and wrap without page overflow.
    page.locator('#tab-sequence-input').focus()
    page.keyboard.press('ArrowRight')
    expect(page.locator('#tab-query-options')).to_be_focused()
    expect(page.locator('#query-options')).to_be_visible()
    page.keyboard.press('End')
    expect(page.locator('#tab-general-options')).to_be_focused()
    page.keyboard.press('Home')
    expect(page.locator('#tab-sequence-input')).to_be_focused()

    # Upload and file drop use local text, validate before replacing, and reset.
    with page.expect_file_chooser() as chooser:
        page.locator('#query-upload-button').click()
    chooser.value.set_files({'name': 'query.fa', 'mimeType': 'text/plain', 'buffer': b'>uploaded\nGGGAGGAGGG\n'})
    expect(page.locator('#query')).to_have_value('>uploaded\nGGGAGGAGGG\n')
    page.locator('#target').evaluate("""node => {
      const transfer = new DataTransfer();
      transfer.items.add(new File(['>dropped\nCCCUCCUCCC\n'], 'target.fa', {type:'text/plain'}));
      node.dispatchEvent(new DragEvent('drop', {dataTransfer:transfer, bubbles:true, cancelable:true}));
    }""")
    expect(page.locator('#target')).to_have_value('>dropped\nCCCUCCUCCC\n')
    page.locator('#query-upload').set_input_files({'name': 'bad.fa', 'mimeType': 'text/plain', 'buffer': b'not a sequence!'})
    expect(page.locator('#input-error')).to_be_visible()
    expect(page.locator('#query')).to_have_value('>uploaded\nGGGAGGAGGG\n')
    # A file read completing after Reset cannot restore the previous input.
    page.evaluate("""() => { window.originalFileText = File.prototype.text; File.prototype.text = function() { return new Promise(resolve => { window.finishRead = () => resolve('>late\nGGGAGGAGGG\n'); }); }; }""")
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

    # Reference prediction, selectable CSV output and independent table visibility.
    fill('outNumber', 10)
    page.locator('#outCsvSort').select_option('E')
    predict()
    expect(rows).to_have_count(10)
    assert rows.first.locator('td').all_text_contents()[1:7] == ['target', '1', '10', 'query', '1', '10']
    assert rows.first.locator('td').all_text_contents()[-2:] == ['-21.40', '0.00']
    page.get_by_role('button', name='Sort by Energy (kcal/mol)', exact=False).click()
    assert float(rows.first.locator('td').nth(9).inner_text()) == -14.3
    page.get_by_role('button', name='Sort by Energy (kcal/mol)', exact=False).click()
    page.locator('#search').fill('-21.4')
    expect(rows).to_have_count(1)
    page.locator('#search').fill('no match')
    expect(page.locator('#no-matches')).to_be_visible()
    page.locator('#search').fill('')
    page.locator('.column-picker summary').click()
    page.locator('#visible-columns input[value="start1"]').uncheck()
    assert 'Target start' not in page.locator('thead').inner_text()
    page.locator('#visible-columns input[value="start1"]').check()
    with page.expect_download() as downloaded:
        page.locator('#download').click()
    csv = Path(downloaded.value.path()).read_text()
    assert csv.splitlines()[0] == 'id1;start1;end1;id2;start2;end2;subseqDP;hybridDP;E'
    assert 'target;1;10;query;1;10;' in csv
    assert 'deltaE' not in csv
    assert not visualization_urls  # Third-party content loads only on selection.
    rows.first.get_by_role('button', name='Visualize').click()
    expect(page.locator('#varri-frame')).to_be_visible()
    expect(page.frame_locator('#varri-frame').get_by_role('heading', name='vaRRI viewer')).to_be_visible()
    expect(rows.first.locator('button')).to_have_attribute('aria-pressed', 'true')
    query = parse_qs(urlparse(page.locator('#varri-frame').get_attribute('src')).query)
    assert query['sequence'] == ['CCCUCCUCCC&GGGAGGAGGG']
    assert query['structure'] == ['((((((((((&))))))))))']
    assert query['seqName1'] == ['target'] and query['seqName2'] == ['query']
    assert query['startIndex1'] == ['1']
    assert query['showRenderingOnly'] == ['true']
    assert 'showRenderingOnly' not in page.locator('#varri-link').get_attribute('href')
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
    expect(page.locator('#results-view')).not_to_be_visible()
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
    fill('query', '>q1\nGGGAGGAGGG\n>q2\nGGGAGGAGGG')
    fill('qIdxPos0', -3)
    fill('tIdxPos0', -5)
    open_group('output-options')
    page.locator('#output-options input[value="qMinE"]').check()
    predict()
    expect(rows).to_have_count(2)
    assert rows.first.locator('td').all_text_contents()[1:7] == ['target', '-5', '5', 'q1', '-3', '7']
    expect(page.locator('#file-list details')).to_have_count(2)
    page.locator('#file-list summary').first.click()
    expect(page.locator('#file-list pre').first).to_contain_text('minE')
    archive = download_archive()
    assert len([name for name in archive if name.startswith('qMinE')]) == 2
    saved = json.loads(archive['parameter.json'])
    assert saved['parameters']['query']['text'].count('>') == 2
    assert saved['parameters']['qIdxPos0'] == -3
    rows.nth(1).get_by_role('button', name='Visualize').click()
    query = parse_qs(urlparse(page.locator('#varri-frame').get_attribute('src')).query)
    assert query['seqName2'] == ['q2']
    assert query['startIndex1'] == ['-5'] and query['startIndex2'] == ['-3']

    # FASTA identity follows native first-token naming, including collisions.
    fill('query', '>q first\nGGGAGGAGGG\n>q second\nGGGAGGAGGG')
    compute.click()
    expect(page.locator('#input-error')).to_contain_text('unique')
    expect(page.locator('#varri-frame')).not_to_be_visible()
    assert page.locator('#varri-frame').get_attribute('src') is None
    expect(page.locator('#download-zip')).not_to_be_visible()

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
    page.locator('#output-options input[value="ED1"]').check()
    predict()
    expect(status).to_contain_text('Evaluation complete')
    expect(rows).to_have_count(2)
    expect(page.locator('thead')).to_contain_text('ED1')
    expect(page.locator('#file-list summary')).to_contain_text('intarna-detailed.txt')
    page.locator('#file-list summary').click()
    expect(page.locator('#file-list pre')).to_contain_text('-21.4')
    with page.expect_download() as report:
        page.locator('#file-list button').click()
    assert '-21.4' in Path(report.value.path()).read_text()
    archive = download_archive()
    assert '-21.4' in archive['intarna-detailed.txt']
    saved = json.loads(archive['parameter.json'])
    assert saved['parameters']['personality'] == 'IntaRNAeval'
    assert saved['parameters']['rri'] == '1||||||||||&1||||||||||:1|||||||||&2|||||||||'
    assert 'seedBP' not in saved['parameters']
    assert 'detailed' in saved['gui']['reports']

    # Visualization stays available when structure/identity columns are deselected.
    open_group('output-options')
    page.get_by_role('button', name='Clear selection', exact=True).click()
    page.locator('#output-options input[value="E"]').check()
    predict()
    rows.first.get_by_role('button', name='Visualize').click()
    expect(page.locator('#varri-frame')).to_be_visible()
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
    page.route('**/wasm/intarna.wasm', lambda route: held_routes.append(route))
    compute.click()
    expect(page.locator('#loading')).to_be_visible()
    page.locator('#abort').click()
    expect(status).to_contain_text('aborted')
    for route in held_routes:
        route.abort()
    page.unroute('**/wasm/intarna.wasm')
    fill('query', '>long_query\n' + 'G' * 1000)
    fill('target', '>long_target\n' + 'C' * 1000)
    compute.click()
    expect(page.locator('#loading-status')).to_contain_text('Predicting RNA', timeout=120000)
    page.keyboard.press('Escape')
    expect(status).to_contain_text('aborted')
    load_example()
    predict()

    fill('query', '>no_pairs\nNNNNNNNNNN')
    compute.click()
    expect(status).to_contain_text('No favorable interactions', timeout=120000)
    expect(rows).to_have_count(0)
    expect(page.locator('#varri-frame')).not_to_be_visible()
    assert download_archive()['intarna-interactions.csv'].count('\n') == 1
    fill('query', '><img/src=x/onerror=alert(1)>\nGGGAGGAGGG')
    predict()
    assert page.locator('tbody img').count() == 0
    assert '<img' in rows.first.inner_text()

    load_example()
    page.route('**/wasm/intarna.js', lambda route: route.abort())
    compute.click()
    expect(status).to_contain_text('failed', timeout=30000)
    expect(page.locator('#loading')).not_to_be_visible()
    expect(compute).to_be_enabled()
    page.unroute('**/wasm/intarna.js')
    predict()
    assert errors == [], errors
    browser.close()
print(f'{args.browser}: real WASM, groups/defaults/versions, prediction/evaluation, FASTA/coordinates, tabs/keyboard, visualization URLs/selection, ZIP replay contents, columns/reports, mobile and abort/recovery passed')

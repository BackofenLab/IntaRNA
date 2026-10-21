"""Exercise the staged browser interface with real WASM at a nested URL."""
import argparse
from pathlib import Path
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
    page.goto(args.url)
    page.wait_for_load_state('networkidle')
    compute = page.locator('#compute')
    status = page.locator('#result-status')
    rows = page.locator('tbody tr')

    def open_group(key):
        group = page.locator('#' + key)
        if not group.evaluate('node => node.open'):
            group.locator('summary').first.click()
        expect(group).to_have_attribute('open', '')
        assert page.locator('.parameter-group[open]').count() == 1

    def fill(key, value):
        parent = page.locator('#' + key).evaluate('node => node.closest(".parameter-group")?.id')
        if parent:
            open_group(parent)
        page.locator('#' + key).fill(str(value))

    def predict():
        compute.click()
        expect(status).to_contain_text('complete:', timeout=120000)
        expect(page.locator('#input-error')).not_to_be_visible()

    expect(page.locator('#version-info')).to_contain_text('ViennaRNA 2.7.0')
    for logo in page.locator('.site-header img').all():
        assert logo.evaluate('node => node.complete && node.naturalWidth > 0')
    for key in ['query-options', 'target-options', 'seed-options', 'helix-options', 'shape-options', 'interaction-options', 'output-options', 'general-options', 'sequence-input']:
        open_group(key)

    # JavaScript validation works even when the offending control is folded.
    fill('query', '>bad\nACGU-ACGU')
    compute.click()
    expect(page.locator('#input-error')).to_contain_text('use only')
    expect(page.locator('#loading')).not_to_be_visible()
    page.locator('#example').click()
    fill('temperature', '101')
    open_group('sequence-input')
    compute.click()
    expect(page.locator('#input-error')).to_contain_text('Temperature')
    page.locator('#reset').click()

    # Reference prediction, selectable CSV output and independent table visibility.
    fill('outNumber', 10)
    page.locator('#outCsvSort').select_option('E')
    predict()
    expect(rows).to_have_count(10)
    assert rows.first.locator('td').all_text_contents()[:6] == ['target', '1', '10', 'query', '1', '10']
    assert rows.first.locator('td').all_text_contents()[-2:] == ['-21.40', '0.00']
    page.get_by_role('button', name='Sort by Energy (kcal/mol)', exact=False).click()
    assert float(rows.first.locator('td').nth(8).inner_text()) == -14.3
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
    if args.screenshot:
        open_group('sequence-input')
        page.screenshot(path=args.screenshot, full_page=True)
    page.set_viewport_size({'width': 390, 'height': 844})
    assert page.evaluate('document.documentElement.scrollWidth <= innerWidth')
    page.set_viewport_size({'width': 1440, 'height': 1100})

    # Multi-FASTA prediction, signed inclusive output and additional files.
    page.locator('#example').click()
    fill('query', '>q1\nGGGAGGAGGG\n>q2\nGGGAGGAGGG')
    fill('qIdxPos0', -3)
    fill('tIdxPos0', -5)
    open_group('output-options')
    page.locator('#output-options input[value="qMinE"]').check()
    predict()
    expect(rows).to_have_count(2)
    assert rows.first.locator('td').all_text_contents()[:6] == ['target', '-5', '5', 'q1', '-3', '7']
    expect(page.locator('#file-list details')).to_have_count(2)
    page.locator('#file-list summary').first.click()
    expect(page.locator('#file-list pre').first).to_contain_text('minE')

    # FASTA identity follows native first-token naming, including collisions.
    fill('query', '>q first\nGGGAGGAGGG\n>q second\nGGGAGGAGGG')
    compute.click()
    expect(page.locator('#input-error')).to_contain_text('unique')

    # Evaluation resets native defaults and removes prediction-only controls.
    fill('query', '>q1\nGGGAGGAGGG\n>q2\nGGGAGGAGGG')
    page.locator('#task').select_option('evaluation')
    expect(page.locator('#seed-options')).not_to_be_visible()
    expect(page.locator('#outNumber')).to_be_disabled()
    fill('rri', '1||||||||||&1||||||||||')
    compute.click()
    expect(page.locator('#input-error')).to_contain_text('exactly one')
    page.locator('#example').click()
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
    fill('rri', '1|||&1||')
    compute.click()
    expect(page.locator('#input-error')).to_contain_text('same nonzero')

    # Startup and synchronous C++ cancellation both allow a fresh prediction.
    page.locator('#task').select_option('prediction')
    page.locator('#example').click()
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
    page.locator('#example').click()
    predict()

    fill('query', '>no_pairs\nNNNNNNNNNN')
    compute.click()
    expect(status).to_contain_text('No favorable interactions', timeout=120000)
    expect(rows).to_have_count(0)
    fill('query', '><img/src=x/onerror=alert(1)>\nGGGAGGAGGG')
    predict()
    assert page.locator('tbody img').count() == 0
    assert '<img' in rows.first.inner_text()

    page.locator('#example').click()
    page.route('**/wasm/intarna.js', lambda route: route.abort())
    compute.click()
    expect(status).to_contain_text('failed', timeout=30000)
    expect(page.locator('#loading')).not_to_be_visible()
    expect(compute).to_be_enabled()
    page.unroute('**/wasm/intarna.js')
    predict()
    assert errors == [], errors
    browser.close()
print(f'{args.browser}: real WASM, groups/defaults/versions, prediction/evaluation, FASTA/coordinates, columns/reports, mobile and abort/recovery passed')

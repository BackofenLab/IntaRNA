"""Browser integration tests against staged assets, including the real WASM engine.

Start an HTTP server serving build-wasm/site, then pass its URL. Optional
--browser selects Firefox; --screenshot saves the populated desktop workspace.
"""
import argparse
from pathlib import Path
from playwright.sync_api import sync_playwright, expect

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("url")
parser.add_argument("--browser", choices=["chromium", "firefox"], default="chromium")
parser.add_argument("--screenshot")
args = parser.parse_args()

with sync_playwright() as playwright:
    browser = getattr(playwright, args.browser).launch(headless=True)
    page = browser.new_page(viewport={"width": 1440, "height": 1100})
    errors = []
    page.on("pageerror", lambda error: errors.append(str(error)))
    page.goto(args.url)
    page.wait_for_load_state("networkidle")
    compute = page.get_by_role("button", name="Compute interactions", exact=True)
    status = page.locator("#result-status")

    # Invalid bases and multiple records must fail before engine startup.
    page.locator("#query").fill(">bad\nACGU-ACGU")
    compute.click()
    expect(page.locator("#input-error")).to_contain_text("use only")
    expect(page.locator("#loading")).not_to_be_visible()
    page.locator("#query").fill(">a\nACGUACGU\n>b\nACGUACGU")
    compute.click()
    expect(page.locator("#input-error")).to_contain_text("exactly one")
    page.get_by_role("button", name="Load example").click()

    # Real RNAlib prediction: orientation, native-reference energy and suboptimals.
    compute.click()
    expect(status).to_contain_text("Prediction complete", timeout=120000)
    rows = page.locator("tbody tr")
    expect(rows).to_have_count(10)
    assert rows.first.locator("td").all_text_contents() == ["query", "1", "10", "target", "1", "10", "-21.40", "0.00"]
    assert page.locator(".table-scroll").evaluate("node => node.scrollWidth") > 0
    for key in ["Query", "Query start", "Query end", "Target", "Target start", "Target end", "Energy (kcal/mol)", "ΔE (kcal/mol)"]:
        page.get_by_role("button", name=f"Sort by {key}", exact=False).first.click()
    # Select energy, then check numeric descending and ascending order.
    page.get_by_role("button", name="Sort by Energy (kcal/mol)", exact=False).click()
    page.get_by_role("button", name="Sort by Energy (kcal/mol)", exact=False).click()
    expect(page.locator('th[aria-sort="descending"]')).to_contain_text("Energy")
    assert float(rows.first.locator("td").nth(6).inner_text()) == -14.3
    page.get_by_role("button", name="Sort by Energy (kcal/mol)", exact=False).click()
    assert float(rows.first.locator("td").nth(6).inner_text()) == -21.4
    page.locator("#search").fill("-21.4")
    expect(rows).to_have_count(1)
    page.locator("#search").fill("no match")
    expect(page.locator("#no-matches")).to_be_visible()
    page.locator("#search").fill("")
    with page.expect_download() as downloaded:
        page.get_by_role("button", name="Download CSV").click()
    csv = Path(downloaded.value.path()).read_text()
    assert "target;1;10;query;1;10;-21.4" in csv
    if args.screenshot:
        page.screenshot(path=args.screenshot, full_page=True)

    # Mobile content stays inside the viewport; the data table scrolls itself.
    page.set_viewport_size({"width": 390, "height": 844})
    assert page.evaluate("document.documentElement.scrollWidth <= innerWidth")
    page.set_viewport_size({"width": 1440, "height": 1100})

    # Deliberately hold WASM download so Abort deterministically covers loading.
    held_routes = []
    page.route("**/wasm/intarna.wasm", lambda route: held_routes.append(route))
    compute.click()
    expect(page.locator("#loading")).to_be_visible()
    page.get_by_role("button", name="Abort computation").click()
    expect(status).to_contain_text("aborted")
    expect(compute).to_be_enabled()
    # Settle any intercepted request after terminating the worker.
    for route in held_routes:
        route.abort()
    page.unroute("**/wasm/intarna.wasm")

    # Wait until actual synchronous C++ starts, then Escape must still work.
    page.locator("#query").fill(">long_query\n" + "G" * 1000)
    page.locator("#target").fill(">long_target\n" + "C" * 1000)
    compute.click()
    expect(page.locator("#loading-status")).to_contain_text("Predicting RNA", timeout=120000)
    page.keyboard.press("Escape")
    expect(status).to_contain_text("aborted")
    page.get_by_role("button", name="Load example").click()
    compute.click()
    expect(status).to_contain_text("Prediction complete", timeout=120000)
    expect(rows).to_have_count(10)

    # Ambiguous bases produce a valid empty prediction, not a fabricated result.
    page.locator("#query").fill(">no_pairs\nNNNNNNNNNN")
    compute.click()
    expect(status).to_contain_text("No favorable interactions", timeout=120000)
    expect(page.locator("#results-view")).not_to_be_visible()

    # Engine-load failures settle the dialog and allow the user to try again.
    page.get_by_role("button", name="Load example").click()
    page.route("**/wasm/intarna.js", lambda route: route.abort())
    compute.click()
    expect(status).to_contain_text("failed", timeout=30000)
    expect(page.locator("#loading")).not_to_be_visible()
    expect(compute).to_be_enabled()
    page.unroute("**/wasm/intarna.js")
    compute.click()
    expect(status).to_contain_text("Prediction complete", timeout=120000)
    assert errors == [], errors
    browser.close()
print(f"{args.browser}: real WASM, validation, sorting/filtering/export, mobile, empty/error and abort/restart passed")

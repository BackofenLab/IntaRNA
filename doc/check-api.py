#!/usr/bin/env python3
"""Check the generated API entry page, its links, and navigation assets."""

import sys
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit


class Page(HTMLParser):
    def __init__(self, path):
        super().__init__()
        self.links = []
        self.anchors = set()
        self.feed(path.read_text(encoding="utf-8"))

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        for key in ("href", "src"):
            if attrs.get(key):
                self.links.append(attrs[key])
        for key in ("id", "name"):
            if attrs.get(key):
                self.anchors.add(attrs[key])


def check(directory):
    for name in ("index.html", "annotated.html", "hierarchy.html",
                 "namespaceIntaRNA.html", "search/searchdata.js", "intarna.css"):
        path = directory / name
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty API output: {path}")

    index = directory / "index.html"
    page = Page(index)
    pages = {index: page}
    for link in page.links:
        url = urlsplit(link)
        if url.scheme or url.netloc:
            continue
        target = directory / unquote(url.path) if url.path else index
        if not target.is_file():
            raise ValueError(f"Broken API entry-page link: {link}")
        if url.fragment and target.suffix == ".html":
            if target not in pages:
                pages[target] = Page(target)
            if unquote(url.fragment) not in pages[target].anchors:
                raise ValueError(f"Broken API entry-page anchor: {link}")

    print(f"API entry page: {len(page.links)} links and assets checked")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: python3 doc/check-api.py <generated-html-directory>")
    try:
        check(Path(sys.argv[1]))
    except (OSError, ValueError) as error:
        sys.exit(str(error))

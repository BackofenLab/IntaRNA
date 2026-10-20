"""Apply narrow dlib 19.24 compatibility fixes in the downloaded ViennaRNA tree.

Modern libc++ only defines char_traits for standard character types. dlib now
uses char32_t upstream. Its old serialization overloads and dependent template
call also need the matching C++ fixes. RNA models and recurrences are untouched.
"""
import sys
from pathlib import Path

root = Path(sys.argv[1]) / "src/dlib-19.24/dlib"
patches = [
    ("unicode/unicode.h", "typedef uint32 unichar;", "typedef char32_t unichar;"),
    ("serialize.h", "    USE_DEFAULT_INT_SERIALIZATION_FOR(wchar_t)",
     "    USE_DEFAULT_INT_SERIALIZATION_FOR(wchar_t)\n    USE_DEFAULT_INT_SERIALIZATION_FOR(char32_t)"),
    ("global_optimization/find_max_global.h", "::template go(std::forward<T>(f),a)",
     "::template go<T>(std::forward<T>(f),a)"),
]
for filename, before, after in patches:
    path = root / filename
    source = path.read_text()
    if after in source:
        continue
    if source.count(before) != 1:
        raise SystemExit(f"Unexpected dlib source: {path}; review compatibility patch")
    path.write_text(source.replace(before, after))

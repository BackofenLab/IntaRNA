# Companion modules and optional executable builder.
SUBDIRS = copomus intarnapvalue bin

dist_bin_SCRIPTS = \
	CopomuS.py \
	IntaRNA-fuse.py \
	IntaRNA_CSV_p-value.py \
	IntaRNA_plotRegions.py

# Standalone scripts import this module from their installed directory.
dist_pythonhelpers_DATA = intarna_tools.py
pythonhelpersdir = $(bindir)
noinst_PYTHON = setup.py

if HAVE_PYTHON
if INSTALL_INTARNAPVALUE
# Install modules using the selected interpreter without fetching dependencies.
# Users install NumPy/SciPy in that interpreter's environment as documented.
install-exec-hook:
	$(PYTHON) -m pip install --ignore-installed --no-deps --no-build-isolation \
		--prefix="$(DESTDIR)$(prefix)" "$(abs_srcdir)"
endif
endif

EXTRA_DIST = README.md DEVELOPMENT.md requirements-output-tools.txt pyproject.toml MANIFEST.in \
	plotRegions.example.png examples/regions.csv examples/regions.params \
	tests/test_output_tools.py tests/test_copomus.py tests/test_intarnapvalue.py \
	tests/test_common_tools.py tests/test_native_tools.py

.PHONY: check-output-tools check-python-tools
check-output-tools: check-python-tools
check-python-tools:
	$(PYTHON) -m unittest discover -s $(srcdir)/tests -v

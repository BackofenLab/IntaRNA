
sourcefiles = \
		__init__.py \
		candidate_filters.py \
		candidate_selectors.py \
		indexing.py \
		measures.py \
		mutation.py \
		mutation_generators.py
		
# install python modules
if HAVE_PYTHON

# copomus installation
if INSTALL_INTARNAPVALUE

copomus_PYTHON = $(sourcefiles)
		
copomusdir = $(pythondir)/copomus

else
noinst_PYTHON = $(sourcefiles)
endif

else

# ensure files are part of distribution if python not available
EXTRA_DIST = $(sourcefiles)

endif

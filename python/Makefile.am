
# sub directories to check for Makefiles
SUBDIRS = copomus intarnapvalue bin


# scripts to be installed in $(bindir) directly

dist_bin_SCRIPTS = \
	CopomuS.py

# auxiliary scripts to be part of distribution 
noinst_PYTHON = setup.py

# install python modules
if HAVE_PYTHON

if INSTALL_INTARNAPVALUE

# ensure pip is installed
# create python directory if non-existing with current prefix
install-exec-hook:
	@ $(PYTHON) -m ensurepip --default-pip
	@ $(MKDIR_P) -p ${pythondir}
	@ export PYTHONPATH=$PYTHONPATH:${pythondir}; \
		$(PYTHON) setup.py install --prefix=${prefix}

endif

endif


# sub directories to check for Makefiles
SUBDIRS = copomus intarnapvalue bin


# install python modules
if HAVE_PYTHON

noinst_PYTHON = setup.py

if INSTALL_INTARNAPVALUE

# ensure pip is installed
# create python directory if non-existing with current prefix
install-exec-hook:
	@ $(PYTHON) -m ensurepip --default-pip
	@ $(MKDIR_P) -p ${pythondir}
	@ export PYTHONPATH=$PYTHONPATH:${pythondir}; \
		$(PYTHON) setup.py install --prefix=${prefix}

endif

# copomus installation
if INSTALL_INTARNAPVALUE

python_PYTHON = \
		CopomuS.py

endif

else

# ensure files are part of distribution if python not available
EXTRA_DIST = *.py

endif

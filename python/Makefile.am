
# sub directories to check for Makefiles
SUBDIRS = intarnapvalue bin

noinst_PYTHON = setup.py


# install python modules
if HAVE_PYTHON
if INSTALL_INTARNAPVALUE

install-exec-hook:
	# ensure pip is installed
	@ $(PYTHON) -m ensurepip --default-pip
	# create python directory if non-existing with current prefix
	@ $(MKDIR_P) -p ${pythondir}
	# compile and install module
	@ export PYTHONPATH=$PYTHONPATH:${pythondir}; \
		$(PYTHON) setup.py install --prefix=${prefix}

endif
endif

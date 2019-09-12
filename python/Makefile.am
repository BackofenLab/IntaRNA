
# sub directories to check for Makefiles
SUBDIRS = intarnapvalue bin

noinst_PYTHON = setup.py


# install python modules
if HAVE_PYTHON
if INSTALL_INTARNAPVALUE
install-exec-hook:
	@ $(PYTHON) -m ensurepip --default-pip
	@ export PYTHONPATH=$PYTHONPATH:${pythondir}; \
		$(PYTHON) setup.py install --prefix=${prefix}
endif
endif

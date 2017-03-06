
# sub directories to check for Makefiles
SUBDIRS = src perl tests .

################################################################
# tests
################################################################

# global test directive # ensure all is compiled (need object files for tests)
.PHONY: tests test
tests: all
	@echo
	@echo " => compiling and running tests ... will take a while ..."
	@echo
	@$(MAKE) check -s -C tests

test: tests

################################################################
# pkg-config
################################################################

if enable_pkg_config
	pkgconfigdir = $(libdir)/pkgconfig
	pkgconfig_DATA = IntaRNA.pc
endif

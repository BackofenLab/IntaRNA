
# sub directories to check for Makefiles
SUBDIRS = src perl .

# global test directive # ensure all is compiled (need object files for tests)
tests: all
	@$(MAKE) tests -s -C tests

test: tests



# sub directories to check for Makefiles
SUBDIRS = src perl tests .

# global test directive # ensure all is compiled (need object files for tests)
tests: all
	@echo
	@echo " => compiling and running tests ... will take a while ..."
	@echo
	@$(MAKE) tests -s -C tests

test: tests


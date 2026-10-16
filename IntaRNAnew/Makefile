CXX ?= g++
AR ?= ar
CPPFLAGS ?= -Iinclude
CXXFLAGS ?= -std=c++23 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Wconversion -Wshadow
DEPFLAGS ?= -MMD -MP
LDFLAGS ?=
LDLIBS ?=
THREAD_FLAGS ?= -pthread

PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin
LIBDIR ?= $(PREFIX)/lib
INCLUDEDIR ?= $(PREFIX)/include
DOCDIR ?= $(PREFIX)/share/doc/IntaRNAnew

LIB_SOURCES := $(wildcard src/*.cpp)
LIB_OBJECTS := $(LIB_SOURCES:.cpp=.o)
TOOLS_LIB_SOURCES := $(wildcard src/tools/*.cpp)
TOOLS_LIB_OBJECTS := $(TOOLS_LIB_SOURCES:.cpp=.o)
APP_OBJECT := app/main.o
TEST_OBJECT := tests/integration.o
FOCUSED_TEST_OBJECT := tests/cli_io_output.o
CONFIG_TEST_OBJECT := tests/config_registry.o
ACCESSIBILITY_TEST_OBJECT := tests/accessibility_oracle.o
FOLDING_TEST_OBJECT := tests/folding_oracle.o
PREDICTOR_TEST_OBJECT := tests/predictor_oracle.o
EXECUTION_TEST_OBJECT := tests/execution_regions.o
MODEL_TEST_OBJECT := tests/model_semantics.o
COMPRESSION_TEST_OBJECT := tests/compression_io.o
OUTPUT_EXECUTION_TEST_OBJECT := tests/output_execution.o
RUNNER_TEST_OBJECT := tests/runner_contracts.o
TOOLS_TEST_OBJECT := tests/tools_contracts.o
PVALUE_EXEC_TEST_OBJECT := tests/pvalue_executable.o
TOOLS_APP_OBJECTS := tools/intarnanew-csv.o tools/intarnanew-stats.o \
	tools/intarnanew-svg.o tools/intarnanew-mutate.o tools/intarnanew-pvalue.o
CORE_LIBRARY := libintarnanew.a
TOOLS_LIBRARY := libintarnanew-tools.a
TOOLS_APPS := intarnanew-csv intarnanew-stats intarnanew-svg intarnanew-mutate \
	intarnanew-pvalue
ALL_OBJECTS := $(LIB_OBJECTS) $(TOOLS_LIB_OBJECTS) $(APP_OBJECT) \
	$(TEST_OBJECT) $(FOCUSED_TEST_OBJECT) $(CONFIG_TEST_OBJECT) \
	$(ACCESSIBILITY_TEST_OBJECT) $(FOLDING_TEST_OBJECT) \
	$(PREDICTOR_TEST_OBJECT) $(EXECUTION_TEST_OBJECT) $(MODEL_TEST_OBJECT) \
	$(COMPRESSION_TEST_OBJECT) $(OUTPUT_EXECUTION_TEST_OBJECT) \
	$(RUNNER_TEST_OBJECT) $(TOOLS_TEST_OBJECT) $(PVALUE_EXEC_TEST_OBJECT) \
	$(TOOLS_APP_OBJECTS)
DEPENDENCIES := $(ALL_OBJECTS:.o=.d)
MISSING_DEPENDENCY_OBJECTS := $(patsubst %.d,%.o,\
	$(filter-out $(wildcard $(DEPENDENCIES)),$(DEPENDENCIES)))
LEGACY_CASES := $(abspath tests/fixtures/parameters)
CONSUMER_SMOKE_SOURCE := tests/consumer/main.cpp

.PHONY: all clean check debug sanitize benchmark legacy-check \
	pvalue-executable-check install install-check force-missing-dependency-rebuild

# Existing objects without compiler dependency files are unsafe after a public
# header change. This can occur when upgrading an older build tree. Rebuild such
# objects once; the next invocation sees their generated .d files and returns to
# normal timestamp-based incremental builds.
ifneq ($(strip $(MISSING_DEPENDENCY_OBJECTS)),)
$(MISSING_DEPENDENCY_OBJECTS): force-missing-dependency-rebuild
endif
force-missing-dependency-rebuild:

all: IntaRNAnew $(TOOLS_APPS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(DEPFLAGS) $(THREAD_FLAGS) -c -o $@ $<

IntaRNAnewBenchmark: benchmarks/compare.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

benchmark: IntaRNAnewBenchmark

$(CORE_LIBRARY): $(LIB_OBJECTS)
	$(RM) $@
	$(AR) rcs $@ $^

$(TOOLS_LIBRARY): $(TOOLS_LIB_OBJECTS)
	$(RM) $@
	$(AR) rcs $@ $^

IntaRNAnew: $(APP_OBJECT) $(CORE_LIBRARY)
	$(CXX) $(LDFLAGS) -o $@ $(APP_OBJECT) $(CORE_LIBRARY) $(LDLIBS) $(THREAD_FLAGS)

IntaRNAnewIntegrationTests: $(TEST_OBJECT) $(CORE_LIBRARY)
	$(CXX) $(LDFLAGS) -o $@ $(TEST_OBJECT) $(CORE_LIBRARY) $(LDLIBS) $(THREAD_FLAGS)

IntaRNAnewCliIoOutputTests: $(FOCUSED_TEST_OBJECT) $(CORE_LIBRARY)
	$(CXX) $(LDFLAGS) -o $@ $(FOCUSED_TEST_OBJECT) $(CORE_LIBRARY) $(LDLIBS) $(THREAD_FLAGS)

IntaRNAnewConfigRegistryTests: $(CONFIG_TEST_OBJECT) $(CORE_LIBRARY)
	$(CXX) $(LDFLAGS) -o $@ $(CONFIG_TEST_OBJECT) $(CORE_LIBRARY) $(LDLIBS) $(THREAD_FLAGS)

IntaRNAnewAccessibilityOracleTests: $(ACCESSIBILITY_TEST_OBJECT) $(CORE_LIBRARY)
	$(CXX) $(LDFLAGS) -o $@ $(ACCESSIBILITY_TEST_OBJECT) $(CORE_LIBRARY) $(LDLIBS) $(THREAD_FLAGS)

IntaRNAnewFoldingOracleTests: $(FOLDING_TEST_OBJECT) $(CORE_LIBRARY)
	$(CXX) $(LDFLAGS) -o $@ $(FOLDING_TEST_OBJECT) $(CORE_LIBRARY) $(LDLIBS) $(THREAD_FLAGS)

IntaRNAnewPredictorOracleTests: $(PREDICTOR_TEST_OBJECT) $(CORE_LIBRARY)
	$(CXX) $(LDFLAGS) -o $@ $(PREDICTOR_TEST_OBJECT) $(CORE_LIBRARY) $(LDLIBS) $(THREAD_FLAGS)

IntaRNAnewExecutionRegionsTests: $(EXECUTION_TEST_OBJECT) $(CORE_LIBRARY)
	$(CXX) $(LDFLAGS) -o $@ $(EXECUTION_TEST_OBJECT) $(CORE_LIBRARY) $(LDLIBS) $(THREAD_FLAGS)

IntaRNAnewModelSemanticsTests: $(MODEL_TEST_OBJECT) $(CORE_LIBRARY)
	$(CXX) $(LDFLAGS) -o $@ $(MODEL_TEST_OBJECT) $(CORE_LIBRARY) $(LDLIBS) $(THREAD_FLAGS)

IntaRNAnewCompressionIoTests: $(COMPRESSION_TEST_OBJECT) $(CORE_LIBRARY)
	$(CXX) $(LDFLAGS) -o $@ $(COMPRESSION_TEST_OBJECT) $(CORE_LIBRARY) $(LDLIBS) $(THREAD_FLAGS)

IntaRNAnewOutputExecutionTests: $(OUTPUT_EXECUTION_TEST_OBJECT) $(CORE_LIBRARY)
	$(CXX) $(LDFLAGS) -o $@ $(OUTPUT_EXECUTION_TEST_OBJECT) $(CORE_LIBRARY) $(LDLIBS) $(THREAD_FLAGS)

IntaRNAnewRunnerContractTests: $(RUNNER_TEST_OBJECT) $(CORE_LIBRARY)
	$(CXX) $(LDFLAGS) -o $@ $(RUNNER_TEST_OBJECT) $(CORE_LIBRARY) $(LDLIBS) $(THREAD_FLAGS)

IntaRNAnewToolsContractTests: $(TOOLS_TEST_OBJECT) $(TOOLS_LIBRARY) $(CORE_LIBRARY)
	$(CXX) $(LDFLAGS) -o $@ $(TOOLS_TEST_OBJECT) $(TOOLS_LIBRARY) $(CORE_LIBRARY) \
		$(LDLIBS) $(THREAD_FLAGS)

intarnanew-csv: tools/intarnanew-csv.o $(TOOLS_LIBRARY) $(CORE_LIBRARY)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) $(THREAD_FLAGS)

intarnanew-stats: tools/intarnanew-stats.o $(TOOLS_LIBRARY) $(CORE_LIBRARY)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) $(THREAD_FLAGS)

intarnanew-svg: tools/intarnanew-svg.o $(TOOLS_LIBRARY) $(CORE_LIBRARY)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) $(THREAD_FLAGS)

intarnanew-mutate: tools/intarnanew-mutate.o $(TOOLS_LIBRARY) $(CORE_LIBRARY)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) $(THREAD_FLAGS)

intarnanew-pvalue: tools/intarnanew-pvalue.o $(TOOLS_LIBRARY) $(CORE_LIBRARY)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) $(THREAD_FLAGS)

IntaRNAnewPvalueExecutableTests: $(PVALUE_EXEC_TEST_OBJECT)
	$(CXX) $(LDFLAGS) -o $@ $< $(LDLIBS)

ifeq ($(OS),Windows_NT)
pvalue-executable-check:
	@echo "SKIP: the p-value executable contract harness requires POSIX process APIs"
else
pvalue-executable-check: intarnanew-pvalue IntaRNAnewPvalueExecutableTests
	./IntaRNAnewPvalueExecutableTests "$(abspath intarnanew-pvalue)"
endif

check: all IntaRNAnewIntegrationTests IntaRNAnewCliIoOutputTests \
	IntaRNAnewConfigRegistryTests \
	IntaRNAnewAccessibilityOracleTests IntaRNAnewFoldingOracleTests \
	IntaRNAnewPredictorOracleTests \
	IntaRNAnewExecutionRegionsTests IntaRNAnewModelSemanticsTests \
	IntaRNAnewCompressionIoTests \
	IntaRNAnewOutputExecutionTests \
	IntaRNAnewRunnerContractTests \
	IntaRNAnewToolsContractTests
	./IntaRNAnewIntegrationTests
	./IntaRNAnewCliIoOutputTests
	./IntaRNAnewConfigRegistryTests
	./IntaRNAnewAccessibilityOracleTests
	./IntaRNAnewFoldingOracleTests
	./IntaRNAnewPredictorOracleTests
	./IntaRNAnewExecutionRegionsTests
	./IntaRNAnewModelSemanticsTests
	./IntaRNAnewCompressionIoTests
	./IntaRNAnewOutputExecutionTests
	./IntaRNAnewRunnerContractTests
	./IntaRNAnewToolsContractTests
	$(MAKE) --no-print-directory pvalue-executable-check
	$(MAKE) --no-print-directory legacy-check

ifeq ($(OS),Windows_NT)
legacy-check:
	@echo "SKIP: the legacy black-box verifier requires POSIX process APIs"
else
legacy-check: IntaRNAnew IntaRNAnewBenchmark
	./IntaRNAnewBenchmark --new "$(abspath IntaRNAnew)" \
		--cases "$(LEGACY_CASES)" --verify-only --expect-cases 17 \
		--skip-without-legacy
endif

debug:
	$(MAKE) clean
	$(MAKE) CXXFLAGS="-std=c++23 -O0 -g3 -Wall -Wextra -Wpedantic -Wconversion -Wshadow" all

sanitize:
	$(MAKE) clean
	$(MAKE) \
		CXXFLAGS="-std=c++23 -O1 -g3 -Wall -Wextra -Wpedantic -Wconversion -Wshadow -fsanitize=address,undefined -fno-omit-frame-pointer" \
		LDFLAGS="-fsanitize=address,undefined" check

install: all $(CORE_LIBRARY) $(TOOLS_LIBRARY)
	install -d "$(DESTDIR)$(BINDIR)" "$(DESTDIR)$(LIBDIR)" \
		"$(DESTDIR)$(INCLUDEDIR)/intarnanew/tools" "$(DESTDIR)$(DOCDIR)/tools"
	install -m 755 IntaRNAnew $(TOOLS_APPS) "$(DESTDIR)$(BINDIR)"
	install -m 644 $(CORE_LIBRARY) $(TOOLS_LIBRARY) "$(DESTDIR)$(LIBDIR)"
	install -m 644 include/intarnanew/*.hpp "$(DESTDIR)$(INCLUDEDIR)/intarnanew"
	install -m 644 include/intarnanew/tools/*.hpp \
		"$(DESTDIR)$(INCLUDEDIR)/intarnanew/tools"
	install -m 644 README.md MATHEMATICAL_AUDIT.md LICENSE \
		"$(DESTDIR)$(DOCDIR)"
	install -m 644 tools/README.md "$(DESTDIR)$(DOCDIR)/tools"

install-check: all
	@set -eu; \
	stage=$$(mktemp -d /tmp/intarnanew-install.XXXXXX); \
	trap 'rm -rf -- "$$stage"' EXIT HUP INT TERM; \
	$(MAKE_COMMAND) --no-print-directory install DESTDIR="$$stage" PREFIX=/usr \
		BINDIR=/usr/bin LIBDIR=/usr/lib INCLUDEDIR=/usr/include \
		DOCDIR=/usr/share/doc/IntaRNAnew; \
	"$$stage/usr/bin/IntaRNAnew" --version >/dev/null; \
	for tool in $(TOOLS_APPS); do "$$stage/usr/bin/$$tool" --help >/dev/null; done; \
	$(CXX) $(CXXFLAGS) -I"$$stage/usr/include" \
		$(CONSUMER_SMOKE_SOURCE) -L"$$stage/usr/lib" \
		-lintarnanew-tools -lintarnanew $(LDLIBS) $(THREAD_FLAGS) \
		-o "$$stage/intarnanew-consumer"; \
	"$$stage/intarnanew-consumer"

clean:
	$(RM) $(ALL_OBJECTS) $(DEPENDENCIES) \
		$(CORE_LIBRARY) $(TOOLS_LIBRARY) IntaRNAnew $(TOOLS_APPS) \
		IntaRNAnewBenchmark \
		IntaRNAnewIntegrationTests IntaRNAnewCliIoOutputTests \
		IntaRNAnewConfigRegistryTests \
		IntaRNAnewAccessibilityOracleTests IntaRNAnewFoldingOracleTests \
		IntaRNAnewPredictorOracleTests \
		IntaRNAnewExecutionRegionsTests IntaRNAnewModelSemanticsTests \
		IntaRNAnewCompressionIoTests \
		IntaRNAnewOutputExecutionTests \
		IntaRNAnewRunnerContractTests \
		IntaRNAnewPvalueExecutableTests \
		IntaRNAnewToolsContractTests

-include $(DEPENDENCIES)

# Logging evaluation: Easylogging++ and ViennaRNA

Evaluation for [issue #273](https://github.com/BackofenLab/IntaRNA/issues/273),
2026-10-16. IntaRNA source baseline:
[`efe5258`](https://github.com/BackofenLab/IntaRNA/tree/efe52586624f0b499c9e040f246bbd5a11d3979b).
This report proposes no change to runtime logging, CLI options or dependencies.

## Recommendation

**Keep the current implementation for now. ViennaRNA is a plausible future
backend, but its built-in logger is not a feature-equivalent replacement.**
Its existing dependency, callback API and logging lock are good reasons to
consider consolidation. Preserving IntaRNA's behavior would still require an
IntaRNA-owned C++ adapter for verbosity, scoped timing, message construction,
formatting and destination selection. Removing Easylogging++ would therefore
trade a bundled general-purpose library for a smaller application-specific
logging layer, rather than remove logging maintenance altogether.

The current ViennaRNA minimum already supports the necessary callback API.
The main decision is whether removing the bundled dependency and simplifying
library initialization justify that adapter and the public API transition.
Styling alone does not justify a replacement: both systems support colors,
but their output and controls differ substantially.

## Scope and evidence

The comparison uses the **bundled Easylogging++ v9.96.7**, including IntaRNA's
local modification, rather than the current upstream Easylogging++ release.
Contrary to the issue's starting description, this checkout does not use it
as a purely header-only dependency:

- [src/Makefile.am](../../src/Makefile.am) compiles `easylogging++.cc` into
  `libeasylogging.a` and installs both that archive and `easylogging++.h`.
- [general.h](../../src/IntaRNA/general.h) includes it from public IntaRNA
  headers; [the executable](../../src/bin/IntaRNA.cpp) and
  [API test runner](../../tests/runApiTests.cpp) initialize its global storage.
- [configure.ac](../../configure.ac) enables performance tracking, suppresses
  creation of a default log file, and selects console coloring. The vendor
  header/source contain 7,683 lines together; this is a source-maintenance
  observation, not a measured build-time or executable-size saving.

ViennaRNA's public [logging documentation](https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/utils/logs.html)
provides the API overview. Implementation-dependent conclusions below refer
to these release-pinned primary sources:

| Release | Sources | Relevance |
| --- | --- | --- |
| 2.7.0 | [log.h](https://github.com/ViennaRNA/ViennaRNA/blob/219394580aec203a9d6f0d5450021e22642d5a83/src/ViennaRNA/utils/log.h), [log.c](https://github.com/ViennaRNA/ViennaRNA/blob/219394580aec203a9d6f0d5450021e22642d5a83/src/ViennaRNA/utils/log.c), [release notes](https://github.com/ViennaRNA/ViennaRNA/blob/219394580aec203a9d6f0d5450021e22642d5a83/CHANGELOG.md) | Minimum supported release; introduced the flexible logger |
| 2.7.2 | [log.h](https://github.com/ViennaRNA/ViennaRNA/blob/1ffec79f5e258896160f7362ced8263450f371dc/src/ViennaRNA/utils/log.h), [log.c](https://github.com/ViennaRNA/ViennaRNA/blob/1ffec79f5e258896160f7362ced8263450f371dc/src/ViennaRNA/utils/log.c) | Additional tested release; `log.c` is identical to 2.7.0 |

IntaRNA now requires **ViennaRNA >= 2.7.0**, including a configure link check
for `vrna_log_lock_set()`. The former >= 2.4.14 range is historical; those older
releases must not be assumed to provide the new logging API. No further
minimum-version increase would be needed for the features evaluated here.
See the existing [ViennaRNA synchronization audit](vienna-accessibility-thread-safety.md)
for the version change and threading contract. This evaluation does not claim
that an unexamined future RNAlib implementation behaves identically.

## What IntaRNA actually uses

An inventory of uncommented logging macro sites in `src/IntaRNA/*.{h,cpp}` and
`src/bin/*.{h,cpp}` at the baseline gives:

| Facility | Sites | Purpose |
| --- | ---: | --- |
| `LOG(ERROR)` | 59 | Argument/input validation and diagnostic failures |
| `LOG(INFO)` | 45 | Ignored/reset options and input interpretation |
| `LOG(WARNING)` | 19 | Invalid/skipped seeds, short sequences, overflows, exceptions |
| `LOG(DEBUG)` | 1 | Unexpected ViennaRNA callback data, inside an IntaRNA debug guard |
| `VLOG(1)` | 11 | Progress, parameters and accessible regions |
| `VLOG(2)` | 30 | Predictor/accessibility details and seed/helix counts |
| `TIMED_FUNC_IF(..., VLOG_IS_ON(9))` | 24 | Scoped timings, emitted through the performance logger |

These are source sites, including templates, not execution counts. There are
no active application `LOG(FATAL)`/`LOG(TRACE)` or performance-checkpoint
sites. For representative examples see [CommandLineParsing](../../src/bin/CommandLineParsing.cpp),
[PredictorMfe2d::predict()](../../src/IntaRNA/PredictorMfe2d.cpp),
[SeedHandlerNoBulge](../../src/IntaRNA/SeedHandlerNoBulge.cpp) and
[AccessibilityVrna](../../src/IntaRNA/AccessibilityVrna.cpp).

The executable configures all registered loggers with `# %level : %msg`.
IntaRNA's own messages go to **stdout**, including errors; prediction output
has its own handlers. `--default-log-file=PATH` sends INFO, WARNING, VERBOSE and
DEBUG messages to the file, while ERROR messages remain on stdout and are
excluded from that file. No default log file is created. Errors are not
silenced by `--default-log-file=/dev/null`. This behavior is documented in the
[README](../../README.md) and implemented before IntaRNA parses its arguments
in [main()](../../src/bin/IntaRNA.cpp).

`-v` and `--verbose` set Easylogging++'s verbosity to **9**, thereby enabling
both progress/detail messages and scoped timings. This is separate from the
ordinary severity levels: INFO and WARNING already appear without `-v`.
The vendor parser also understands numeric verbosity/module arguments, but
IntaRNA's Boost options register only `verbose,v` and `default-log-file` for
logging, and reject unregistered options. Vendor features such as `-vmodule`
are therefore not established IntaRNA CLI requirements. The C++ logging API
can still select verbosity programmatically.

The timings produce INFO messages such as `Executed [function] in [duration]`.
They measure elapsed work with Easylogging++'s scope object, not merely the
wall-clock timestamp of a message. Replacing them with RNAlib's timestamp
option would lose this information.

## Capability comparison

| Requirement | Current IntaRNA/Easylogging++ | ViennaRNA 2.7.0 / 2.7.2 | Consequence of replacement |
| --- | --- | --- | --- |
| Severity | INFO, WARNING, ERROR and guarded DEBUG used | DEBUG, INFO, WARNING, ERROR, CRITICAL, plus SILENT threshold; default is ERROR | Straightforward mapping, but leaving defaults would hide current INFO/WARNING messages |
| Verbosity | Separate VLOG scale; sites at 1/2 and timing gate at 9 | No verbose scale or verbosity field in log events | Keep an IntaRNA gate/category; mapping all verbosity to DEBUG would merge it with diagnostics |
| Scoped duration | 24 conditional scope timers; named performance logger | Optional timestamp and file/line metadata, no duration tracker | Implement a scoped timer and preserve its gate/format |
| C++ message construction | Stream insertion supports strings, numbers and project types such as `IndexRangeList` | `printf`-style variadic C API | Convert call sites or retain a lazy C++ front end; guard disabled messages before constructing them |
| Format | Configurable format, currently `# LEVEL : message` | Built-in padded `[LEVEL] message`; optional time and source location | Exact output requires a custom callback; `vrna_log_options_set()` cannot set a format template |
| Destination | Per-level file/console configuration; ERROR exception described above | One default `FILE *`; callbacks have independent minimum levels | A callback must implement the exact routing; simply setting one file pointer sends errors there too |
| File ownership | Vendor opens/manages the path and flushes logs | Caller supplies/manages `FILE *`; built-in sink flushes each message | Adapter owns opening, failure handling, lifetime and closing |
| CLI | Existing options plus initialization before argument parsing | RNAlib API has no IntaRNA option parsing | Retain IntaRNA flags; ViennaRNA tools' `--log-*` options are not inherited by linking RNAlib |
| Other RNAlib diagnostics | A separate logger with separate destinations/thresholds | Application messages and RNAlib messages share global state | Unified formatting/routing becomes possible, but changes can also reveal or suppress dependency messages |
| Threading | Vendor has optional locking; this build does not define `ELPP_THREAD_SAFE`; some CLI/library sites use OpenMP critical sections | Lock hook covers a whole log dispatch, including callbacks | Useful consolidation, but initialize before any application logging and obey callback restrictions below |
| Embedding | Public headers expose Easylogging++; consumers initialize it | Public C API has process-global settings and callbacks | Simpler consumer startup is possible, but global RNAlib configuration needs an explicit ownership policy |

Both systems have more generic features than IntaRNA currently needs. Named
loggers, arbitrary configuration files, rotation, asynchronous logging and
rate-limited logging should not be counted as mandatory gaps simply because
the vendor supports them. The `performance` logger is used; arbitrary named
application channels are not. `LogDetailedCrashReason` is set in `main()`,
but the normal configure flags do not enable `ELPP_FEATURE_CRASH_LOG` or
`ELPP_FEATURE_ALL`; a full crash/backtrace facility is not an active feature
that this migration must reproduce.

## Styling and stream behavior

The implementations are in the bundled [easylogging++.cc](../../src/easylogging++.cc)
(`LogBuilder::convertToColoredOutput()`, `OS::termSupportsColor()`,
`DefaultLogDispatchCallback::dispatch()`) and ViennaRNA's pinned `log.c`
(`log_default()`, `get_log_level_color()`).

| Detail | IntaRNA/Easylogging++ | ViennaRNA built-in sink |
| --- | --- | --- |
| Plain text examples | `# INFO : message`, `# WARNING : message`, `# VERBOSE : message` | `[INFO]    message`, `[WARNING] message`; no VERBOSE label |
| Colored region | Entire line, including message | Severity label/padding; optional source location has its own gray color |
| Colors | INFO cyan, WARNING yellow, ERROR/FATAL red, DEBUG green; VERBOSE unchanged | INFO bold blue, WARNING bold yellow, ERROR bold red, DEBUG bold cyan, CRITICAL/FATAL bold magenta |
| When enabled | IntaRNA defaults to coloring; `--disable-log-coloring` removes the flag | Depends on how RNAlib was compiled (`VRNA_WITHOUT_TTY_COLORS`) and whether its destination is a TTY |
| Terminal detection | Vendor checks a fixed `TERM` whitelist, not whether stdout is a TTY | `isatty(fileno(destination))` |
| Redirection | stdout can still contain ANSI escapes in a pipe when `TERM=xterm`; the explicit log-file sink stays plain | Redirected/file output stays plain with the tested builds |
| Runtime control | Application controls the color flag | No public per-call/runtime color switch in the audited logging API; use a custom sink for application-owned control |

A move to RNAlib's default styling would visibly change labels, spacing,
colors and stdout/stderr behavior. Its TTY detection is an advantage for
redirected output, but that improvement can also be made independently of a
backend migration. IntaRNA's `--disable-log-coloring` cannot control an already
compiled RNAlib's built-in coloring. A callback can preserve that IntaRNA
build option without requiring users to rebuild ViennaRNA.

## Risks that an adapter must address

1. **Global state and initialization.** `VrnaHandler` currently installs only
   the RNAlib lock, once, preserving streams, levels and user callbacks.
   Reconfiguring the global sink/threshold in each handler would change that
   contract for embedding applications. IntaRNA logs during startup before a
   handler exists; an RNAlib backend would need earlier lock setup. Configure
   callbacks and destinations before workers begin, and remove only callbacks
   owned by IntaRNA after workers stop. `vrna_log_reset()` removes callbacks
   and the installed lock; it must not be used as routine application setup.
2. **Callback locking and lifetime.** RNAlib dispatch invokes callbacks under
   its logging lock. A callback must not log back into RNAlib, take the setup
   mutex while another path holds it and logs, or let C++ exceptions escape.
   Keep callback data and output streams alive for the full logging lifetime.
   A future bridge from RNAlib to Easylogging++ would inherit these constraints
   and would not by itself repair Easylogging++'s native-thread safety.
3. **Callback filtering and variadic arguments.** The default severity
   threshold and `QUIET` option do not suppress registered callbacks, which
   have their own thresholds. Turning off the default sink is not a global
   mute. In both audited `log.c` files, `log_default()` consumes `event->params`
   with `vfprintf()` before `log_v()` passes that same event to callbacks;
   the dispatcher does not make a fresh `va_copy` per destination. Do not
   assume default output plus consuming callbacks is a safe fan-out design.
   For an IntaRNA-owned sink, quiet the default output and format from a local
   `va_copy`/`va_end` without consuming the shared list. Coexistence with
   third-party callbacks still needs separate validation and an ownership
   policy; no upstream fix or multi-callback safety claim is made here.
4. **Severity is not termination.** `vrna_log_error()` and
   `vrna_log_critical()` return after logging. The older
   `vrna_message_error()` exits the process in both audited implementations.
   It is not a replacement for IntaRNA's nonterminating `LOG(ERROR)` calls
   followed by normal validation/exception handling. IntaRNA disables abort
   on Easylogging++ fatal logs, although it has no active FATAL sites.
5. **Source and build compatibility.** Removing the vendor affects installed
   headers/archive, consumer initialization, tests, linking instructions and
   `configure.ac`'s source-directory marker, not just the executable. Update
   those together in a future migration. A logging facade must retain lazy
   verbosity checks and safe formatting of project types; C varargs alone do
   not provide stream insertion's type support. RNAlib's debug macro is also
   subject to `VRNA_LOG_NO_DEBUG_RNALIB`, independently of IntaRNA's debug guard.

The current OpenMP critical sections do not establish a general native-thread
logging guarantee, particularly for scoped timer destruction and other
unguarded paths. A synchronized replacement could improve that boundary;
it would still need tests, not just removal of the old critical sections.

## Options and decision criteria

| Option | Benefits | Costs / limitations | Assessment |
| --- | --- | --- | --- |
| Retain Easylogging++ | Preserves present CLI, formatting, timing and consumer behavior | Two logging systems and vendored code remain; existing threading/color limitations remain | Appropriate outcome for this evaluation-only issue |
| Replace calls directly with RNAlib macros | Removes vendor, reuses an existing dependency and lock | Loses verbose categories/timers; changes format, destination, defaults and C++ construction | Do not pursue as a behavior-preserving change |
| IntaRNA facade with an RNAlib callback sink | Can remove vendor and consumer initialization; can unify RNAlib diagnostics and centralize locking | IntaRNA must own verbosity, timers, formatting, routing, lifecycle and embedding policy | Viable future project if those maintenance benefits justify the work |
| Bridge RNAlib diagnostics into existing Easylogging++ | Can make dependency diagnostics look consistent with application messages | Does not remove vendor; introduces callback/lock/global-state concerns | Consider only if unified diagnostics are themselves a requested feature |

No runtime speedup, build-time reduction or binary-size reduction has been
measured. The reasons to consider replacement are dependency maintenance,
consumer ergonomics and a unified diagnostic policy. Acceptance of a future
migration should require preserving default and verbose output, scoped timing,
file/error routing, early parse failures, no-color builds, pipe/TTY behavior,
parallel logging and embedding with existing RNAlib callbacks. Any deliberate
change to those behaviors should be documented as such.

## Validation and reproducible RNAlib probe

The analysis used GCC 14.4.0 with `-std=c++23` on Linux. Isolated probes compiled
against the bundled Easylogging++ source and ViennaRNA **2.7.0 and 2.7.2**
confirmed:

- Easylogging++ with the executable's configuration writes to stdout; `-v`
  enables both VLOG(1/2) messages and a scoped timing INFO message. With a
  configured log file, those messages go to the file while ERROR stays on
  stdout. Files have no ANSI escapes; captured stdout still has them with
  `TERM=xterm` and the color flag enabled.
- Both RNAlib releases suppress INFO/WARNING at the default ERROR threshold,
  print to stderr, and return after a critical log. At INFO threshold they
  print INFO/WARNING too. Their quiet mode still invokes the separately
  registered INFO callback. The legacy error function exits with status 1.
- Pipe captures of RNAlib output have no ANSI escapes. A pseudo-terminal
  capture with `TERM=xterm` has the severity colors listed above. `log.c` is
  byte-identical between the examined release source trees.

This is a focused API/style check, not a logging throughput benchmark,
concurrency stress test, macOS/Windows validation or a full IntaRNA rebuild.
No application source, build flags or logging implementation changed.

Save the following as `logging-probe.cpp`. It intentionally uses only one
callback with the default sink quiet and runs in its own process; its reset
must not be copied into an embedding application's initialization path.

```cpp
#include <cassert>
#include <cstdarg>
#include <cstdio>
#include <cstring>
extern "C" {
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/utils/log.h>
}
static void capture(vrna_log_event_t *event, void *data) {
  ++*static_cast<int *>(data);
  va_list args;
  va_copy(args, event->params);
  std::vprintf(event->format_string, args);
  va_end(args);
  std::putchar('\n');
}
int main(int argc, char **argv) {
  vrna_log_reset();
  assert(vrna_log_level() == VRNA_LOG_LEVEL_ERROR);
  if (argc > 1 && !std::strcmp(argv[1], "legacy-error")) {
    vrna_message_error("legacy error");
    return 99;
  }
  if (argc > 1 && !std::strcmp(argv[1], "callback")) {
    int calls = 0;
    vrna_log_options_set(VRNA_LOG_OPTION_QUIET);
    vrna_log_cb_add(capture, &calls, nullptr, VRNA_LOG_LEVEL_INFO);
    vrna_log_debug("hidden debug");
    vrna_log_info("info=%d", 7);
    vrna_log_warning("warning=%d", 8);
    assert(calls == 2);
    const auto removed = vrna_log_cb_remove(capture, &calls);
    assert(removed == 1);
  } else {
    vrna_log_info("hidden info");
    vrna_log_warning("hidden warning");
    vrna_log_error("visible error");
    vrna_log_level_set(VRNA_LOG_LEVEL_INFO);
    vrna_log_info("visible info");
    vrna_log_warning("visible warning");
    vrna_log_critical("critical returns");
  }
  std::puts("returned");
}
```

With the desired ViennaRNA installation selected by `PKG_CONFIG_PATH` and a
C++23 compiler selected by `CXX`, run. For a custom dependency prefix, also
make its link/runtime libraries discoverable (for example via `LIBRARY_PATH`
and `LD_LIBRARY_PATH` on Linux):

```bash
"${CXX:-c++}" -std=c++23 logging-probe.cpp \
  $(pkg-config --cflags --libs RNAlib2) -o logging-probe
./logging-probe >default.stdout 2>default.stderr
./logging-probe callback >callback.stdout 2>callback.stderr
./logging-probe legacy-error >legacy.stdout 2>legacy.stderr
# The last command intentionally exits with status 1.
```

The first mode exits 0, writes `returned` to stdout and four visible severity
messages to stderr. The callback mode exits 0, writes `info=7`, `warning=8`,
`returned` to stdout, and leaves stderr empty. The legacy mode writes only its
error to stderr and exits 1. Run the first mode in a terminal to inspect color;
redirecting its output should remove the ANSI sequences in the tested builds.

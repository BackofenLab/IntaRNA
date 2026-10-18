# Historical compiler probes and compatibility experiment

These files were recovered from the local issue #244 validation directory.
They accompany the [historical C++23 assessment](../5-cxx23-justification.md)
of commit `0a4568ad6e3da52f935221ff0da98a4a171e9833`.

- `cxx23-library-probe.cpp` exercises the library facilities in that revision's
  configure probe. It checks availability, rather than application performance.
- `cxx17-library-probe.cpp` checks the retained `try_emplace` requirement.
- `set-aside-cxx17.patch` and its notes preserve an abandoned compatibility
  proposal against the historical source. They are not applied by the build
  or CI and have not been rebased or validated against current master.
  Statements in the saved notes describe the proposal, not the supported
  compiler policy. The recorded issue decision retained C++23.

These were local investigation artifacts rather than proposed production
changes. The full audit was left untracked after its shorter conclusion was
posted on issue #244. The compatibility patch was explicitly set aside;
its filename and the saved consensus record establish that it was not part
of the final decision. They are committed now to avoid losing that context
when the local workspace is removed.

For a compiler supporting the corresponding language mode:

```sh
c++ -std=c++23 -pedantic-errors cxx23-library-probe.cpp -o /tmp/intarna-cxx23-probe
/tmp/intarna-cxx23-probe
c++ -std=c++17 -pedantic-errors cxx17-library-probe.cpp -o /tmp/intarna-cxx17-probe
/tmp/intarna-cxx17-probe
```

The local Git attribute preserves normal unified-diff context whitespace in
the archived patch; it does not alter whitespace checks for production code.

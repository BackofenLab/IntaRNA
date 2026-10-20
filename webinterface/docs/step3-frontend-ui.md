# Step 3: sequence input interface (#286)

## Design

The interface is a small scientific workspace: two equally weighted strand
editors above one results region. Query uses teal (`#087f83`), target uses blue
(`#315fb8`), with navy text (`#17324d`), pale blue background (`#f1f7fa`), white
working surfaces and rust (`#aa4b1d`) for errors. Trebuchet/Segoe UI system fonts
keep the page self-contained; only nucleotide input and structure data use
monospace. The paired-strand mark is the single visual motif. All working text
is left aligned. A generic dashboard with metric cards was rejected because
sequence comparison, input direction and coordinate interpretation matter here.

```text
IntaRNA                                      Documentation / Source
RNA–RNA interaction prediction               Runs on this device
[ Query 5′ → 3′                 ][ Target 5′ → 3′                ]
[ FASTA editor                 ][ FASTA editor                  ]
Maximum results / Load example                          Compute
Results / search / export
[ sortable interaction coordinates and energy                  ]
```

## Input and state

`protocol.js` accepts one FASTA record or one plain sequence per input. Wrapped
lines, lowercase and DNA T are normalized to uppercase RNA U. The only allowed
nucleotides are A/C/G/U/T/N; N does not match a complementary base. Empty names,
multiple/misplaced headers, invalid symbols and lengths outside 7–2,000 nt are
rejected. Query length × target length is limited to 1,000,000 for browser
memory. The result count is an integer from 1 to 100. These limits are also
enforced in the worker and do not change the native CLI.

The example is a short complementary pair with a seven-base seed. Predictions
use the native defaults (heuristic prediction, Turner energy model, accessibility
included, 37 °C, seed length 7) and allow overlapping suboptimal interactions.
All sequences are entered 5′ to 3′; the CLI handles reverse query accessibility.

The loading dialog reports initialization versus calculation without inventing
percentages. Abort and Escape terminate the worker and restore the form. Inputs
are disabled during a prediction; the results are cleared when inputs change or
a new prediction begins. Errors and empty output have explicit messages.
Diagnostics are available separately from CSV. Focus returns to Compute after
the dialog closes. There are visible keyboard focus styles, associated labels,
live status/error regions, reduced-motion support and a single-column mobile
layout. Browser tests cover validation, real computation, abort/restart, worker
load failure, responsive overflow and result interactions.

Serve the staged directory over HTTP; module workers cannot load through
`file://`. No frontend bundler, external font, CDN or API service is needed.

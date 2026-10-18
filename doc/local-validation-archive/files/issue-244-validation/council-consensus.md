# Issue #244 consensus record

Source audited: `0a4568ad6e3da52f935221ff0da98a4a171e9833`.

Three independent reviews covered technical requirements, scientific and
performance evidence, and maintenance/support policy. All recommended
retaining C++23 as the current implementation baseline, while distinguishing
purposeful use from scientific necessity or demonstrated performance gains.

The technical review requested a precision change from “Lowering only the
configure requirement would break these uses” to “Changing only the configure
requirement would not make the current source compatible with earlier
standards.” The final text includes that correction.

All three reviewers explicitly approved the exact final wording saved in
`issue-comment.md`. The root reviewer also approved it. The public comment
does not mention the review process or agents. No source/configuration change
or issue closure is part of the decision.

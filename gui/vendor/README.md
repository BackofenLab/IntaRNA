# Plotly browser bundle

`plotly-cartesian.min.js` is the unmodified browser bundle from
`plotly.js-cartesian-dist-min` 4.1.2 (MIT). It provides the scatter, bar and
heatmap traces used for additional IntaRNA output visualizations.
`plotly.css` is a static copy of its runtime plot/modebar rules. Plotly normally
inserts those rules with CSSOM; the bundle's supported `no-inline-styles` sentinel
suppresses that insertion while the same-origin copy retains the GUI's strict
`style-src 'self'` content-security policy.

- npm integrity: `sha512-3UkyAT7/U8f3sOLavLnUjA6mqFmnqEe246GZhxT0J3TT9oZ0mh6SgwsL1yo/3/4xfmRSDFdGOxe+EsvT2XiE0g==`
- bundle SHA-256: `aaf12ff1d4671f7e93f6df4291651ab4ca6b76c0ea5f1576ec1030e3c2819a20`
- static CSS SHA-256: `3bcf2679b557067d8a3492b74092e946544d1cb8b2a72f278bcba95df59afd33`
- source: <https://www.npmjs.com/package/plotly.js-cartesian-dist-min/v/4.1.2>

Keep the adjacent `plotly-LICENSE.txt` when updating the bundle.

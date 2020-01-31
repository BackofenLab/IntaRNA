
# Auxiliary R scripts of the IntaRNA package

- [IntaRNA_CSV_p-value.R](#intarnacsvpvaluer)
- [IntaRNA_plotRegions.R](#intrarnaplotregionsr)


# `IntaRNA_CSV_p-value.R`
### p-value estimates of IntaRNA energies from present energy distribution

For genome-wide target sRNA predictions, we assume the set of energy values 
predicted for all target sufficiently diverse to be used as a background energy
model for minimum energies of putative target sequences. Thus, we can fit a
generalized extreme value (GEV) distribution to the data that is subsequently
used to estimate p-values for each energy.

The `IntaRNA_CSV_p-value.R` script takes an IntaRNA CSV output file (assumed to be
sufficiently large and sane enough for GEV fitting) to compute respective
p-value estimates. The output consists of the input table extended with a
`p-value` column. If no output file is given or in- and output file names are 
equal, the input file is overwritten! 
You can (optionally) specify the column name for which p-values are to be 
estimated. Example calls are given below.

```bash
# overwriting the input file with p-value-extended table
Rscript --vanilla IntaRNA_CSV_p-value.R IntaRNA-output.csv
# creating a new output file for p-value extension
Rscript --vanilla IntaRNA_CSV_p-value.R IntaRNA-output.csv IntaRNA-output-with-pValue.csv
# computing p-values for normalized energies (has to be present in file IN.csv)
Rscript --vanilla IntaRNA_CSV_p-value.R IN.csv IN-pValue.csv E_norm
```


# `IntaRNA_plotRegions.R`
### Visualization of RRI-covered regions

To visualize sequences' regions covered by RNA-RNA interactions predicted by
IntaRNA, you can use `IntaRNA_plotRegions.R` by providing the following arguments (in 
the given order)

1. CSV-IntaRNA output file (semicolon separated) covering the columns `start,end,id`
  with suffix `1` or `2` to plot target or query regions, respectively
2. `1` or `2` to select whether to plot target or query regions
3. output file name with a file-format-specific suffix from `.pdf`, `.png`, 
  `.svg`, `.eps`, `.ps`, `.jpeg`, `.tiff`

An example is given below, when calling
```bash
Rscript --vanilla IntaRNA_plotRegions.R pred.csv 1 plotRegions.example.png
```

with `pred.csv` containing
```
id1;start1;end1;id2;start2;end2
b0001;266;273;query;116;123
b0002;204;231;query;85;111
b0003;229;262;query;96;125
b0004;265;300;query;10;38
b0005;281;295;query;5;22
```

will produce the output
![IntaRNA_plotRegions.R example](plotRegions.example.png)


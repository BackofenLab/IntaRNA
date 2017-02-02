[![Build Status](https://travis-ci.org/BackofenLab/IntaRNA.svg?branch=master)](https://travis-ci.org/BackofenLab/IntaRNA)

# IntaRNA version 1.*
Efficient target prediction incorporating accessibility of interaction sites

**Motivation**: During the last few years, several new small regulatory RNAs (sRNAs) have been discovered in bacteria. Most of them act as post-transcriptional regulators by base pairing to a target mRNA, causing translational repression or activation, or mRNA degradation. Numerous sRNAs have already been identified, but the number of experimentally verified targets is considerably lower. Consequently, computational target prediction is in great demand. Many existing target prediction programs neglect the accessibility of target sites and the existence of a seed, while other approaches are either specialized to certain types of RNAs or too slow for genome-wide searches.

**Results:** We introduce INTARNA, a new general and fast approach to the prediction of RNA–RNA interactions incorporating accessibility of target sites as well as the existence of a user-definable seed. We successfully applied INTARNA to the prediction of bacterial sRNA targets and determined the exact locations of the interactions with a higher accuracy than competing programs. 

## Contribution

Feel free to contribute to this project by wirting [Issues](https://github.com/BackofenLab/IntaRNA/issues) with feature requests or bug reports.

## Cite
If you use IntaRNA, please cite our [article](http://bioinformatics.oxfordjournals.org/content/24/24/2849):
```
doi: 10.1093/bioinformatics/btn544
```


# IntaRNA [![releases](https://img.shields.io/github/tag/BackofenLab/IntaRNA.svg)](https://github.com/BackofenLab/IntaRNA/releases)  [![Bioconda](https://anaconda.org/bioconda/intarna/badges/version.svg)](https://anaconda.org/bioconda/intarna) [![Docker Repository on Quay](https://quay.io/repository/biocontainers/intarna/status "Docker Repository on Quay")](https://quay.io/repository/biocontainers/intarna) [![Build Status](https://travis-ci.org/BackofenLab/IntaRNA.svg?branch=master)](https://travis-ci.org/BackofenLab/IntaRNA)

**Efficient RNA-RNA interaction prediction incorporating accessibility and 
seeding of interaction sites**

During the last few years, several new small regulatory RNAs 
(sRNAs) have been discovered in bacteria. Most of them act as post-transcriptional 
regulators by base pairing to a target mRNA, causing translational repression 
or activation, or mRNA degradation. Numerous sRNAs have already been identified, 
but the number of experimentally verified targets is considerably lower. 
Consequently, computational target prediction is in great demand. Many existing 
target prediction programs neglect the accessibility of target sites and the 
existence of a seed, while other approaches are either specialized to certain 
types of RNAs or too slow for genome-wide searches.

IntaRNA, developed by
[Prof. Backofen's bioinformatics group at Freiburg University](http://www.bioinf.uni-freiburg.de),
is a general and fast approach to the 
prediction of RNA-RNA interactions incorporating both the accessibility of 
interacting sites 
as well as the existence of a user-definable seed interaction. We successfully applied 
IntaRNA to the prediction of bacterial sRNA targets and determined the exact 
locations of the interactions with a higher accuracy than competing programs. 

For testing or ad hoc use of IntaRNA, you can use its webinterface at the

**==> [Freiburg RNA tools IntaRNA webserver](http://rna.informatik.uni-freiburg.de/IntaRNA/) <==**


## Contribution

Feel free to contribute to this project by writing 
[Issues](https://github.com/BackofenLab/IntaRNA/issues) 
with feature requests, bug reports, or just contact messages.

## Citation

If you use IntaRNA, please cite our articles

- [IntaRNA 2.0: enhanced and customizable prediction of RNA–RNA interactions](http://dx.doi.org/10.1093/nar/gkx279)
  Martin Mann, Patrick R. Wright, and Rolf Backofen, 
  Nucleic Acids Research, 45 (W1), W435–W439, 2017, DOI(10.1093/nar/gkx279).
- [CopraRNA and IntaRNA: predicting small RNA targets, networks and interaction domains](http://dx.doi.org/10.1093/nar/gku359)
  Patrick R. Wright, Jens Georg, Martin Mann, Dragos A. Sorescu, Andreas S. Richter, Steffen Lott, Robert Kleinkauf, Wolfgang R. Hess, and Rolf Backofen,
  Nucleic Acids Research, 42 (W1), W119-W123, 2014, DOI(10.1093/nar/gku359).
- [IntaRNA: efficient prediction of bacterial sRNA targets incorporating target site accessibility and seed regions](http://dx.doi.org/10.1093/bioinformatics/btn544)
  Anke Busch, Andreas S. Richter, and Rolf Backofen, 
  Bioinformatics, 24 no. 24 pp. 2849-56, 2008, DOI(10.1093/bioinformatics/btn544).



<br /><br /><br /><br />
<a name="doc" />

# Documentation

## Overview

The following topics are covered by this documentation:

- [Installation](#install)
  - [IntaRNA via conda](#instconda)
  - [IntaRNA docker container](#instdocker)
  - [Dependencies](#deps)
  - [Cloning from github](#instgithub)
  - [Source code distribution](#instsource)
  - [Microsoft Windows installation](#instwin)
  - [OS X installation with homebrew](#instosx)
- [Usage and Parameters](#usage)
  - [Just run ...](#defaultRun)
  - [Prediction modes, their features and emulated tools](#predModes)
  - [Interaction restrictions](#interConstr)
  - [Seed constraints](#seed)
  - [Explicit seed input](#seedExplicit)
  - [SHAPE reactivity data to enhance accessibility computation](#shape)
  - [Output modes](#outmodes)
  - [Suboptimal RNA-RNA interaction prediction and output restrictions](#subopts)
  - [Energy parameters and temperature](#energy)
  - [Additional output files](#outFiles)
    - [Minimal energy profiles](#profileMinE)
    - [Spot probability profiles](#profileSpotProb) using partition functions
    - [Minimal energy for all intermolecular index pairs](#pairMinE)
    - [Interaction probabilities for interaction spots of interest](#spotProb)
    - [Accessibility and unpaired probabilities](#accessibility)
      - [Local versus global unpaired probabilities](#accLocalGlobal)
      - [Constrain regions to be accessible or blocked](#accConstraints)
      - [Read/write accessibility from/to file or stream](#accFromFile)
  - [Multi-threading and parallelized computation](#multithreading)
- [Library for integration in external tools](#lib)



<br /><br /><br /><br />
<a name="install" />

# Installation

<br /><br />
<a name="instconda" />

## IntaRNA via conda (bioconda channel)

The most easy way to locally install IntaRNA is via conda using the 
[bioconda](https://bioconda.github.io/) 
channel (linux only). This way, you will install a pre-built IntaRNA binary along
with all dependencies.
Follow
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/intarna/README.html)
to get detailed information or run
```bash
conda install -c bioconda intarna
```
if you are using bioconda already.

<br /><br />
<a name="instdocker" />

## IntaRNA docker container (via QUAY)

An [IntaRNA docker container](https://quay.io/repository/biocontainers/intarna) 
([?](https://www.docker.com/)) is provided from the bioconda package via 
[Quay.io](https://quay.io/). This provides
you with an encapsulated IntaRNA installation.


<br /><br />
<a name="deps" />

## Dependencies

If you are going to compile IntaRNA from source, ensure you meet the following
dependencies:

- compiler supporting C++11 standard and OpenMP
- [boost C++ library](http://www.boost.org/) version >= 1.50.0 
  (ensure the following libraries are installed for development (not just runtime libraries!); or install all e.g. in Ubuntu via package `libboost-all-dev`)
    - libboost_regex
    - libboost_program_options
    - libboost_filesystem
    - libboost_system
- [Vienna RNA package](http://www.tbi.univie.ac.at/RNA/) version >= 2.4.8
- if [cloning from github](#instgithub): GNU autotools (automake, autoconf, ..)

Also used by IntaRNA, but already part of the source code distribution (and thus
not needed to be installed separately):

- [Catch](https://github.com/philsquared/Catch) test framework
- [Easylogging++](https://github.com/easylogging/easyloggingpp) logging framework



<br /><br />
<a name="instgithub" />

## Cloning *Source code* from github (or downloading ZIP-file)

The data provided within the github repository
(or within the `Source code` archives provided at the  
[IntaRNA release page](https://github.com/BackofenLab/IntaRNA/releases))
is no complete distribution and
lacks all system specifically generated files. Thus, in order to get started with 
a fresh clone of the IntaRNA source code repository you have to run the GNU autotools 
to generate all needed files for a proper `configure` and `make`. To this end,
we provide the helper script `autotools-init.sh` that can be run as shown in the following.
```bash
# call aclocal, automake, autoconf
bash ./autotools-init.sh
```
Afterwards, you can continue as if you would have downloaded an 
[IntaRNA package distribution](#instsource).

<br /><br />
<a name="instsource" />

## IntaRNA package distribution (e.g. `intaRNA-2.0.0.tar.gz`)

When downloading an IntaRNA package distribution (e.g. `intaRNA-2.0.0.tar.gz`) from the 
[IntaRNA release page](https://github.com/BackofenLab/IntaRNA/releases), you should 
first ensure, that you have all [dependencies](#deps) installed. If so, you can
simply run the following (assuming `bash` shell).
```bash
# generate system specific files (use -h for options)
./configure
# compile IntaRNA from source
make
# run tests to ensure all went fine
make tests
# install (use 'configure --prefix=XXX' to change default install directory)
make install
# (alternatively) install to directory XYZ
make install prefix=XYZ
```

If you installed one of the dependencies in a non-standard directory, you have 
to use the according `configure` options:
- `--with-vrna` : the prefix where the Vienna RNA package is installed
- `--with-boost` : the prefix where the boost library is installed

Note, the latter is for instance the case if your `configure` call returns an
error message as follows:
```[bash]
checking whether the Boost::System library is available... yes
configure: error: Could not find a version of the library!
```
In that case your boost libraries are most likely installed to a non-standard 
directory that you have to specify either using `--with-boost` or just the 
library directory via `--with-boost-libdir`. 

<br /><br />
<a name="instwin" />

## Microsoft Windows installation

### ... from source

IntaRNA can be compiled, installed, and used on a Microsoft Windows system when
e.g. using [Cygwin](https://www.cygwin.com/) as 'linux emulator'. Just install
Cygwin with the following packages:

- *Devel*:
   - make
   - gcc-g++
   - autoconf
   - automake
   - pkg-config
- *Libs*:
   - libboost-devel
- *Perl*:
   - perl

and follow either [install from github](#instgithub) or 
[install from package](#instsource).

### ... using pre-compiled binaries

For some releases, we also provide precompiled binary packages for Microsoft Windows at the
[IntaRNA release page](https://github.com/BackofenLab/IntaRNA/releases) 
that enable 'out-of-the-box' usage. If you
want to use them:
- [download](https://github.com/BackofenLab/IntaRNA/releases) the according ZIP archive and extract
- open a [Windows command prompt](https://www.lifewire.com/how-to-open-command-prompt-2618089)
- [run IntaRNA](#usage) 

*Note*, these binaries come without any waranties, support or what-so-ever!
They are just an offer due to according user requests.

If you do not want to work within the IntaRNA directory or don't want to provide
the full installation path with every IntaRNA call, you should add the installation
directory to your [`Path` System variable](http://www.computerhope.com/issues/ch000549.htm)
(using a semicolon `;` separator).



<br /><br />
<a name="instosx" />

## OS X installation with homebrew (thanks to Lars Barquist)

If you do not want to or can use the pre-compiled binaries for OS X available from 
[bioconda](https://anaconda.org/bioconda/intarna), you can compile `IntaRNA` 
locally.

The following wraps up how to build `IntaRNA-2.0.2` under OS X (Sierra 10.12.4) using homebrew.

First, install homebrew! :)

```[bash]
brew install gcc --without-multilib
```

`--without-multilib` is necessary for OpenMP multithreading support -- note 
OS X default `gcc`/`clang` doesn't support OpenMP, so we need to install standard 
`gcc`/`g++`

```[bash]
brew install boost --cc=gcc-6
```

`--cc=gcc-6` is necessary to build `boost` with standard `gcc`, rather than the 
default bottle which appears to have been built with the system `clang`. 
Brew installs `gcc`/`g++` as `/usr/local/bin/gcc-VERSION` by default to avoid 
clashing with the system's `gcc`/`clang`. `6` is the current version as of 
writing, but may change.

```[bash]
brew install viennarna
brew install doxygen
```

Download and extract the IntaRNA source code package (e.g. `intaRNA-2.0.2.tar.gz`) from the [release page](releases/).

```[bash]
./configure CC=gcc-6 CXX=g++-6
```

This sets up makefiles to use standard `gcc`/`g++` from brew, which will 
need an update to the appropriate compiler version if not still `6`. 
You might also want to
set `--prefix=INSTALLPATH` if you dont want to install IntaRNA globally.


```[bash]
Make
make tests
make install
```




<br /><br /><br /><br />
<a name="usage" />

# Usage and parameters

IntaRNA comes with a vast variety of ways to tune or enhance *YOUR* RNA-RNA prediction.
To this end, different [prediction modes](#predModes) are implemented that allow
to balance predication quality and runtime requirement. Furthermore, it is 
possible to define 
[interaction restrictions](#interConstr),
[seed constraints](#seed), 
[explicit seed information](#seedExplicit), 
[SHAPE reactivity constraints](#shape), 
[output modes](#outmodes),
[suboptimal enumeration](#subopts), 
[energy parameters, temperature](#energy),
and the [accessibility](#accessibility) handling. If you are doing high-throughput
computations, you might also want to consider [multi-threading support](#multithreading).

For ad hoc usage you can use the 
[Freiburg RNA tools IntaRNA webserver](http://rna.informatik.uni-freiburg.de/IntaRNA/)
(with limited parameterization).



<br /><br />
<a name="defaultRun" />

## Just run ...

If you just want to start and are fine with the default parameters set, 
you only have to provide two RNA sequences, 
a (long) target RNA (using `-t` or `--target`) and a (short) query RNA
(via `-q` or `--query`), in 
[IUPAC RNA encoding](#https://en.wikipedia.org/wiki/Nucleic_acid_notation).
You can either directly input the sequences
```bash
# running IntaRNA with direct sequence input
# call : IntaRNA -t CCCCCCCCGGGGGGGGGGGGGG -q CCCCCCC

target
             9     15
             |     |
  5'-CCCCCCCC       GGGGGGG-3'
             GGGGGGG
             |||||||
             CCCCCCC
          3'-       -5'
             |     |
             7     1
query

interaction energy = -10.7116 kcal/mol

```

or provide (multiple) sequence(s) in [FASTA-format](#https://en.wikipedia.org/wiki/FASTA_format).
It is possible to provide either file input or to read the FASTA input from the
STDIN stream.
 
```bash
# running IntaRNA with FASTA files
IntaRNA -t myTargets.fasta -q myQueries.fasta
# reading query FASTA input from stream via pipe
cat myQueries.fasta | IntaRNA -q STDIN -t myTargets.fasta
```

If you are working with large FASTA input files, e.g. covering a whole 
transcriptome, you can restrict the prediction to a subset of the input 
sequences using the `--qSet` or `--tSet` parameter as shown in the following.

```bash
# restrict prediction to the second load of 100 target sequences 
IntaRNA -t myTranscriptome.fasta --tSet=101-200 -q myQuery.fasta
```

Nucleotide encodings different from `ACGUT` are rewritten as `N` and the respective
positions are not considered to form base pairs (and this ignored).
Thymine `T` encodings are replaced by uracil `U`, since a `ACGU`-only 
energy model is used.

For a list of general program argument run `-h` or `--help`. For a complete
list covering also more sophisticated options, run `--fullhelp`.





<br /><br />
<a name="predModes" />

## Prediction modes, their features and emulated tools

For the prediction of *minimum free energy interactions*, the following modes
and according features are supported and can be set via the `--mode` parameter.
The tiem and space complexities are given for the prediction of two sequences
of equal length *n*.

| Features | Heuristic `--mode=H` | Exact-SE `--mode=S` | Exact `--mode=E` |
| -------- | :------------------: | :-----------------: | :--------------: |
| Time complexity (prediction only) | O(*n*^2) | O(*n*^4) | O(*n*^4) |
| Space complexity | O(*n*^2) | O(*n*^2) | O(*n*^4) |
| [Seed constraint](#seed) | x | x | x |
| [Explicit seeds](#seedExplicit) | x | x | x |
| [SHAPE reactivity constraint](#shape) | x | x | x |
| No [seed constraint](#seed) | x | x | x |
| Minimum free energy interaction | not guaranteed | x | x |
| Overlapping [suboptimal interactions](#subopts) | x | x | x |
| Non-overlapping [suboptimal interactions](#subopts) | x | - | x |

Note, due to the low run-time requirement of the heuristic prediction mode
(`--mode=H`), heuristic IntaRNA interaction predictions are widely used to screen
for interaction in a genome-wide scale. If you are more interested in specific
details of an interaction site or of two relatively short RNA molecules, you 
should investigate the exact prediction mode (`--mode=S`, or `--mode=E`
if non-overlapping suboptimal prediction is required).

Given these features, we can emulate and extend a couple of RNA-RNA interaction
tools using IntaRNA.

**TargetScan** and **RNAhybrid** are approaches that predict the interaction hybrid with 
minimal interaction energy without consideration whether or not the interacting 
subsequences are probably involved involved in intramolecular base pairings. Furthermore,
no seed constraint is taken into account.
This prediction result can be emulated (depending on the used prediction mode) 
by running IntaRNA when disabling both the seed constraint
as well as the accessibility integration using
```bash
# prediction results similar to TargetScan/RNAhybrid
IntaRNA [..] --noSeed --qAcc=N --tAcc=N
```
We *add seed-constraint support to TargetScan/RNAhybrid-like computations* by removing the 
`--noSeed` flag from the above call.

**RNAup** was one of the first RNA-RNA interaction prediction approaches that took the 
accessibility of the interacting subsequences into account while not considering the seed feature. 
IntaRNA's exact prediction mode is eventually an alternative implementation when disabling
seed constraint incorporation. Furthermore, the unpaired probabilities used by RNAup to score
the accessibility of subregions are covering the respective overall structural ensemble for each
interacting RNA, such that we have to disable accessibility computation based on local folding (RNAplfold)
using
```bash
# prediction results similar to RNAup
IntaRNA --mode=S --noSeed --qAccW=0 --qAccL=0 --tAccW=0 --tAccL=0
```
We *add seed-constraint support to RNAup-like computations* by removing the 
`--noSeed` flag from the above call.






<br /><br />
<a name="interConstr" />

## Interaction restrictions

The predicted RNA-RNA interactions can be enhanced if additional
knowledge is available. To this end, IntaRNA provides different options to 
restrict predicted interactions.

A general most general restriction is the maximal energy (inversely related to 
stability) an RNA-RNA interaction is allowed to have. Per default, a reported
interaction should have a negative energy (<0) to be energetically favorable.
This report barrier can be altered using `--outMaxE`. For suboptimal interaction
restriction, please refer to [suboptimal interaction prediction](#subopts) section.

If you are only interested in predictions for highly accessible regions, i.e. 
with a high probability to be unpaired, you can use the `--outMinPu` parameter.
If given, each individual position of the interacting subsequences has to have
an unpaired probability reaching at least the given value. This significantly
increases prediction time but will exclude predictions where the formation of
the interaction (intermolecular base pairing) replaces intramolecular base 
pairing (where the latter will cause low unpaired probabilities for the 
respective positions).

Furthermore, the region where interactions are supposed to occur can be restricted
for target and query independently. To this end, a list of according 
subregion-defining index pairs
can be provided using `--qRegion` and `--tRegion`, respectively. The indexing 
starts with 1 and should be in the format `from1-end1,from2-end2,..` using
integers. Note, if you want to have predictions individually for each region
combination (rather than just the best for each query-target combination) you
want to add `--outPerRegion` to the call.

If you are dealing with very long sequences it might be useful to use the
*automatic identification of accessible regions*, which dramatically reduces
runtime and memory consumption of IntaRNA since predictions are only done for
individual regions and not for the whole sequence. Here, we use a 
heuristic approach that finds and ignores subregions that are unlikely to form
an interaction, resulting in a decomposition of the full sequence range into
intervals of accessible regions. It can be enabled by providing the maximal 
length of the resulting intervals via the parameters `--qRegionLenMax` and
`--tRegionLenMax`.<br />
More specifically, starting from the full 
sequence's index range, the algorithm iteratively identifies in every too-long
range the window with highest ED value (penalty for non-accessibility). To
this end, it uses windows of length `--seedBP` to find subsequences where it is
most unlikely that a seed might be formed. This window is removed from the range,
which results in two shorter ranges. If a range is shorter than `--seedBP`, it
is completely removed.

Finally, it is possible to restrict the overall length an interaction is allowed
to have. This can be done independently for the query and target sequence using
`--qIntLenMax` and `--tIntLenMax`, respectively. Both default to the full sequence
length (by setting both to 0).







<br /><br />
<a name="seed" />

## Seed constraints

For different types of RNA-RNA interactions it was found that experimentally
verified interactions were showing a short and compact subinteraction of high
stability (= low energy). It was hypothesized that these regions are among the
first formed parts of the full RNA-RNA interaction, and thus considered as the
*seed* of the overall interaction.

Based on this observation, RNA-RNA interaction predictors were enhanced by
incorporating such seed constraints into their prediction pipeline, i.e. a
reported interaction has to feature at least one seed. Typically, 
a seed is defined as a short subinteraction of 5-8 consecutive base pairs that 
are not enclosing any unpaired nucleotides (or if so only very few).

IntaRNA supports the definition of such seed constraints and adds further
options to even more constrain the seed selection. The list of options is given 
by 

- `--seedBP` : the number of base pairs the seed has to show
- `--seedMaxUP` : the maximal overall number of unpaired bases within the seed
- `--seedQMaxUP` : the maximal number of unpaired bases within the query's seed region
- `--seedTMaxUP` : the maximal number of unpaired bases within the target's seed region
- `--seedMaxE` : the maximal overall energy of the seed (to exclude weak seed interactions)
- `--seedMinPu` : the minimal unpaired probability of each seed region in query and target
- `--seedQRange` : a list of index intervals where a seed in the query is allowed
- `--seedTRange` : a list of index intervals where a seed in the target is allowed

Alternatively, you can set 

- `--seedTQ` : to specify [explicit seed interactions](seedExplicit)

Seed constraint usage can be globally disabled using the `--noSeed` flag.




<br /><br />
<a name="seedExplicit" />

## Explicit seed input

Some experiments provide hints or explicit knowledge about the seed or
even provide details about some intermolecular base pairs formed between two RNAs.
This information can be incorporated into IntaRNA predictions by providing
*explicit seed information*. To this end, the `--seedTQ` parameter can be used.
It takes a comma-separated list of seed string encodings in the format 
`startTbpsT&startQbpsQ`, which is in the same format as the IntaRNA `hybridDB`
output (see below), i.e. e.g. `--seedTQ='4|||.|&7||.||'` 
(ensure you quote the seed encoding to avoid a shell interpretation of the pipe symbol '|') 
to encode a seed interaction like 
the following 
```bash
target
             4    8
             |    |
      5'-AAAC    C UGGUUUGG-3'
             AC C C
             || | |
             UG G G
      3'-GGUU  U   CCCACAAA-5'
             |    |
            11    7
query
```
If several or alternative seeds are known, you can provide all as a 
comma-separated list and IntaRNA will consider all interactions that cover at
least one of them.



<br /><br />
<a name="shape" />

## SHAPE reactivity data to enhance accessibility computation

For some RNA sequences, experimental reactivity data is available that can be
used to guide/help the structure and thus accessibility prediction for the RNA
molecule. IntaRNA supports such data by interfacing the Vienna RNA package
capabilities for SHAPE reactivity data incorporation, see 
Lorenz et al. ([2015](https://doi.org/10.1093/bioinformatics/btv523), [2016](https://dx.doi.org/10.1186%2Fs13015-016-0070-z)) or the
[RNAfold manpage](https://www.tbi.univie.ac.at/RNA/RNAfold.1.html).

The SHAPE reactivity data can be provided via file using `--qShape` or
`--tShape` for query or target sequence, respectively. 
Independently for each, it is possible
to define the methods to be used to convert the data into pseudo energies and
pairing probabilities. The respective IntaRNA arguments are
`--qShapeMethod`|`--tShapeMethod`
and `--qShapeConversion`|`--tShapeConversion`, which mimics the according 
tool arguments in the Vienna RNA package (see e.g. the 
[RNAfold manpage](https://www.tbi.univie.ac.at/RNA/RNAfold.1.html)).




<br /><br />
<a name="outmodes" />

## Output modes

The RNA-RNA interactions predicted by IntaRNA can be provided in different
formats. The style is set via the argument `--outMode` and the different modes
will be discussed below.

Furthermore, it is possible to define *where to output*, i.e. using `--out` 
you can either name a file or one of the stream names `STDOUT`|`STDERR`. Note,
any string not matching one of the two stream names is considered a file name.
The file will be overwritten by IntaRNA!

Besides interaction output, you can set the verbosity of computation information
using the `-v` or `--verbose` arguments. To reduce the output to a minimum, you
can redirect all logging output of user information, warnings or verbose output
to a specific file using `--default-log-file=LOGFILENAME`. 
If you are not interested in any logging output, redirect it to nirvana via 
`--default-log-file=/dev/null`. Note, error output is not redirected and always
given on standard output streams.


<a name="outModeNormal" />

### Standard RNA-RNA interaction output with ASCII chart

The standard output mode `--outMode=D` provides a detailed ASCII chart of the 
interaction together with its overall interaction energy. 
For an example see the [Just run ...](#defaultRun) section.


<a name="outModeDetailed" />

### Detailed RNA-RNA interaction output with ASCII chart

Using `--outMode=D`, a detailed ASCII chart of the interaction together with
various interaction details will be provided. An example is given below.

```bash
# call: IntaRNA -t AAACACCCCCGGUGGUUUGG -q AAACACCCCCGGUGGUUUGG --outMode=D --seedBP=4

target
             5      12
             |      |
      5'-AAAC   C    UGGUUUGG-3'
             ACC CCGG
             ||| ||||
             UGG GGCC
      3'-GGUU   U    CCCACAAA-5'
             |      |
            16      9
query

interaction seq1   = 5 -- 12
interaction seq2   = 9 -- 16

interaction energy = -2.78924 kcal/mol
  = E(init)        = 4.1
  + E(loops)       = -13.9
  + E(dangleLeft)  = -0.458042
  + E(dangleRight) = -0.967473
  + E(endLeft)     = 0.5
  + E(endRight)    = 0
  + ED(seq1)       = 3.91068
  + ED(seq2)       = 4.0256
  + Pu(seq1)       = 0.00175516
  + Pu(seq2)       = 0.00145659

seed seq1   = 9 -- 12
seed seq2   = 9 -- 12
seed energy = -1.4098 kcal/mol
seed ED1    = 2.66437 kcal/mol
seed ED2    = 2.66437 kcal/mol
seed Pu1    = 0.0132596
seed Pu2    = 0.0132596
```
Position annotations start indexing with 1 at the 5'-end of each RNA. 
`ED` values are the energy penalties for reduced [accessibility](#accessibility)
and `Pu` denote unpaired probabilities of the respective interacting subsequences.

<a name="outModeCsv" />

### Customizable CSV RNA-RNA interaction output

IntaRNA provides via `--outMode=C` a flexible interface to generate RNA-RNA 
interaction output in CSV format (using `;` as separator). Note, target sequence
information is listed with index `1` while query sequence information is given
by index `2`.

```bash
# call: IntaRNA -t AAACACCCCCGGUGGUUUGG -q AAACACCCCCGGUGGUUUGG --outMode=C --noSeed --outOverlap=B -n 3
id1;start1;end1;id2;start2;end2;subseqDP;hybridDP;E
target;4;14;query;4;14;CACCCCCGGUG&CACCCCCGGUG;((((...((((&))))...))));-4.14154
target;5;16;query;5;16;ACCCCCGGUGGU&ACCCCCGGUGGU;(((((.((.(((&))))).)).)));-4.04334
target;1;14;query;4;18;AAACACCCCCGGUG&CACCCCCGGUGGUUU;(((((((...((((&))))...)))).)));-2.94305
```
For each prediction, a row in the CSV is generated.

Using the argument `--outCsvCols`, the user can specify what columns are 
printed to the output using a comma-separated list of colIds. Available colIds 
are

- `id1` : id of first sequence (target)
- `id2` : id of second sequence (query)
- `seq1` : full first sequence
- `seq2` : full second sequence
- `subseq1` : interacting subsequence of first sequence
- `subseq2` : interacting subsequence of second sequence
- `subseqDP` : hybrid subsequences compatible with hybridDP
- `subseqDB` : hybrid subsequences compatible with hybridDB
- `start1` : start index of hybrid in seq1
- `end1` : end index of hybrid in seq1
- `start2` : start index of hybrid in seq2
- `end2` : end index of hybrid in seq2
- `hybridDP` : hybrid in VRNA dot-bracket notation
- `hybridDB` : hybrid in dot-bar notation
- `E` : overall hybridization energy
- `ED1` : ED value of seq1
- `ED2` : ED value of seq2
- `Pu1` : probability to be accessible for seq1
- `Pu2` : probability to be accessible for seq2
- `E_init` : initiation energy
- `E_loops` : sum of loop energies (excluding E_init)
- `E_dangleL` : dangling end contribution of base pair (start1,end2)
- `E_dangleR` : dangling end contribution of base pair (end1,start2)
- `E_endL` : penalty of closing base pair (start1,end2)
- `E_endR` : penalty of closing base pair (end1,start2)
- `E_hybrid` : energy of hybridization only = E - ED1 - ED2
- `E_norm` : length normalized energy = E / ln(length(seq1)*length(seq2))
- `E_hybridNorm` : length normalized energy of hybridization only = E_hybrid / ln(length(seq1)*length(seq2))
- `seedStart1` : start index of the seed in seq1
- `seedEnd1` : end index of the seed in seq1
- `seedStart2` : start index of the seed in seq2
- `seedEnd2` : end index of the seed in seq2
- `seedE` : overall hybridization energy of the seed only (excluding rest)
- `seedED1` : ED value of seq1 of the seed only (excluding rest)
- `seedED2` : ED value of seq2 of the seed only (excluding rest)
- `seedPu1` : probability of seed region to be accessible for seq1
- `seedPu2` : probability of seed region to be accessible for seq2

Using `--outCsvCols ''`, all available columns are added to the output.

Energies are provided in unit *kcal/mol*, probabilities in the interval [0,1].
Position annotations start indexing with 1.

The `hybridDP` format is a dot-bracket notation as e.g. generated by **RNAup**.
Here, for each target sequence position within the interaction, 
a '.' represents a position not involved
in the interaction while a '(' marks an interacting position. For the query
sequence this is done analogously but using a ')' for interacting positions.
Both resulting strings are concatenated by a separator '&' to yield a single
string encoding of the interaction's base pairing details.

The `hybridDB` format is similar to the `hybridDP` but also provides site information. 
Here, a bar '|' is used in both base pairing encodings (which makes it a 'dot-bar encoding'). 
Furthermore, each interaction string is prefixed 
with the start position of the respective interaction site.

In the following, an altered CSV output for the example from above is generated.
```bash
# call: IntaRNA --outCsvCols=Pu1,Pu2,subseqDB,hybridDB -t AAACACCCCCGGUGGUUUGG -q AAACACCCCCGGUGGUUUGG --outMode=C --noSeed --outOverlap=B -n 3
Pu1;Pu2;subseqDB;hybridDB
0.00133893;0.00133893;4CACCCCCGGUG&4CACCCCCGGUG;4||||...||||&4||||...||||
0.00134094;0.00134094;5ACCCCCGGUGGU&5ACCCCCGGUGGU;5|||||.||.|||&5|||||.||.|||
0.00133686;0.0013368;1AAACACCCCCGGUG&4CACCCCCGGUGGUUU;1|||||||...||||&4||||...||||.|||
```


<a name="outModeV1" />

### Backward compatible IntaRNA v1.* output

If your scripts/whatever is tuned to the old IntaRNA v1.* output, you can use

- `--outMode=1` : IntaRNA v1.* normal output
- `--outMode=O` : IntaRNA v1.* detailed output (former `-o` option)

Note, for for IntaRNA v1.* output, currently *no multi-threading computation* is available!


<br /><br />
<a name="subopts" />

## Suboptimal RNA-RNA interaction prediction and output restrictions

Besides the identification of the optimal (e.g. minimum-free-energy) RNA-RNA 
interaction, IntaRNA enables the enumeration of suboptimal interactions. To this
end, the argument `-n N` or `--outNumber=N` can be used to generate up to `N`
interactions for each query-target pair (including the optimal one). Note, the
suboptimal enumeration is increasingly sorted by energy.

Note: suboptimal interaction enumeration is not exhaustive! That is, for each
interaction site (defined by the left- and right-most intermolecular base pair)
only the best interaction is reported! In heuristic prediction mode (default
mode of IntaRNA), this is even less exhaustive, since only for each left-most
interaction boundary one interaction is reported!

Furthermore, it is possible to *restrict (sub)optimal enumeration* using

- `--outMaxE` : maximal energy for any interaction reported
- `--outDeltaE` : maximal energy difference of suboptimal interactions' energy
  to the minimum free energy interaction
- `--outOverlap` : defines if an where overlapping of reported interaction sites
  is allowed (Note, IntaRNA v1.* used implicitly the 'T' mode):
  - 'N' : no overlap neither in target nor query allowed for reported interactions
  - 'B' : overlap allowed for interacting subsequences for both target and query
  - 'T' : overlap allowed for interacting subsequences in target only 
  - 'Q' : overlap allowed for interacting subsequences in query only 





<br /><br />
<a name="energy" />

## Energy parameters and temperatures

The selection of the correct temperature and energy parameters is cruicial for
a correct RNA-RNA interaction prediction. To this end, various settings are 
supported by IntaRNA.

The temperature can be set via `--temperature=C`to set a temperature `C` in 
degree Celsius. Note, this is important especially for predictions within plants
etc., since the default temperature is 37°C.

The energy model used can be specified using the `--energy` parameters using

- 'B' for base pair maximization similar to the Nussinov intramolecular structure prediction.
  Here, each base pair contributes an energy term of `-1` independently of its
  structural or sequence context. This mode is mainly useful for study or teaching 
  purpose.
- 'V' enables *Nearest Neighbor Model* energy computation similar to the Zuker
  intramolecular structure prediction using the Vienna RNA package routines. 
  Within this model, the energy contribution of a base
  pair depends on its directly enclosed (neighbored) basepair and the subsequence(s)
  involved. Different energy parameter sets have been experimentally derived
  in the last decades. Since IntaRNA makes use of the energy evaluation routines
  of the Vienna RNA package, all parameter sets from the Vienna RNA package are
  available for RNA-RNA interaction prediction. Per default, the default parameter
  set of the linked Vienna RNA package version is used. You can change the parameter
  set using the `--energyVRNA` parameter as explained below.

If Vienna RNA package is used for energy computation (`--energy=V`), per default
the default parameter set of the linked Vienna RNA package is used (e.g. the
set `Turner04` for VRNA 2.3.0). If you want to use a different parameter set, you
can provide an according parameter file via `--energVRNA=MyParamFile`. The 
following example exemplifies the use of the old `Turner99` parameter set as
used by IntaRNA v1.*.
```bash
# IntaRNA v1.* like energy parameter setup
IntaRNA --energyVRNA=/usr/local/share/Vienna/rna_turner1999.par --seedMaxE=999
```

To increase prediction quality and to reduce the computational complexity, the
number of unpaired bases between intermolecular base pairs is restricted
(similar to internal loop length restrictions in the Zuker algorithm). The
upper bound can be set independent for the query and target sequence via
`--qIntLoopMax` and `--tIntLoopMax`, respectively, and default to 16.






<br /><br />
<a name="outFiles" />

## Additional output files

IntaRNA v2 enables the generation of various additional information in dedicated
files/streams. The generation of such output is guided by an according (repeated)
definition of the `--out` argument in combination with one of the following
argument prefixes (case insensitive) that have to be colon-separated to the
targeted file/stream name:

- `qMinE:`|`tMinE:` the query/target's minimal interaction energy profile (CSV format), respectively
- `pMinE:` the minimal interaction energy for all pairs of query-target index pairs (CSV format)
- `qAcc:`|`tAcc:` the [query/target's ED accessibility values](#accessibility) (RNAplfold-like format), respectively
- `qPu:`|`tPu:` the [query/target's unpaired probabilities](#accessibility) (RNAplfold format), respectively

Note, for *multiple sequences* in FASTA input, the provided file names
are suffixed with with `-s` and the according sequence's number (where indexing
starts with 1) within the FASTA input to generate an individual file for each sequence.
For `pMinE` output, if multiple query-target combinations are to be reported, the
file name is suffixed with `-q#t#` (where `#` denotes the according sequence number
within the input.


<br />
<a name="profileMinE" />

### Minimal energy profiles

To get a more global view of possible interaction sites for a pair of interacting
RNAs, one can generate the *minimal energy profile* for each sequence (independently).

For instance, to generate the target's profile, add the following to your IntaRNA
call: `--out=tMinE:MYPROFILEFILE.csv`. For the query's profile, use `--out=qMinE:..` respectively.
This will produce an according CSV-file (`;` separated) with the according minimal
energy profile data that can be visualized with any program of your liking.

In the following, such an output was visualized using R:
```R
d <- read.table("MYPROFLEFILE.csv", header=T, sep=";");
plot( d[,1], d[,3], xlab="sequence index", ylab="minimal energy", type="l", col="blue", lwd=2)
abline(h=0, col="red", lty=2, lwd=2)
```

![Minimal interaction energy profile of an RNA](/doc/figures/profile-minE.png?raw=true "Minimal interaction energy profile of an RNA")

This plot reveals two less but still stable (*E* below 0) interaction sites beside the
mfe interaction close to the 5'-end of the molecule.



<br />
<a name="profileSpotProb" />

### Spot probability profiles

Similarly to (minimal energy profiles)[#profileMinE], it is also possible to
compute position-wise probabilities how likely a position is covered by an
interaction, i.e. its *spot probability*. To the end, we compute for each
position $i$ the partition function $Zi$ of all interactions covering $i$.
Given the overall partition function *Z* including all possible interactions, 
the position-speficit spot probability for *i* is given by *Zi/Z*.

Such profiles can be generated using `--out=qSpotProb:MYPROFILEFILE.csv` or 
`--out=tSpotProb:...` for the query/target sequence respectively and independently.



<br />
<a name="pairMinE" />

### Minimal energy for all intermolecular index pairs

To investigate how stable RNA-RNA interactions are distributed for a given pair
of RNAs, you can also generate the minimal energy for all intermolecular index
pairs using `--out=pMinE:MYPAIRMINE.csv`. This generates a CSV file (`;`separated)
holding for each index pair the minimal energy of any interaction covering this
index combination or `NA` if no covers it at all.

This information can be visualized with your preferred program. In the following,
the provided R call is used to generate a heatmap visualization of the 
RNA-RNA interaction possibilities.

```R
# read data, skip first column, and replace NA and E>0 values with 0
d <- read.table("MYPAIRMINE.csv",header=T,sep=";");
d <- d[,2:ncol(d)];
d[is.na(d)] = 0;
d[d>0] = 0;
# plot
image( 1:nrow(d), 1:ncol(d), as.matrix(d), col = heat.colors(100), xlab="index in sequence 1", ylab="index in sequence 2");
box();
```

The following plot (for the [minimal energy profile](#profileMinE) example from
above) reveals, that the alternative stable (*E*<0) interactions all involve the
mfe-site in the second sequence and are thus less likely to occure. 

![Minimal interaction energy index pair information](/doc/figures/pair-minE.png?raw=true "Minimal interaction energy index pair information")



<br />
<a name="spotProb" />

### Interaction probabilities for interaction spots of interest

For some research questions, putative regions of interactions are known from
other sources and it is of interest to study the effect of competitive binding
or other scenarios that might influence the accessibility of the interacting
RNAs (e.g. refer to [SHAPE data](#shape) or 
[structure/accessibility constraints](#accConstraints)).

To this end, one can specify the spots of interest by intermolecular index pairs,
e.g. using `5&67` to encode the fifth target RNA position (first number of the
encoding) and the 67th query RNA
position (second number of the encoding). Note, indexing starts with 1. 
Multiple spots can be provided as comma-separated list. The list in
concert with an output stream/file name (colon-separated) can be passed via the
`--out` argument using the `spotProb:` prefix, e.g.

```[bash]
IntaRNA ... --out=spotProb:5&67,33&12:mySpotProbFile.csv
```

The reported probability is the ratio of according partition functions. That is,
for each interaction `I` that respects all input constraints and has an energy 
below 0 (or set `--outMaxE` value) the respective Boltzmann weight `bw(I)` 
is computed by `bw(I) = exp( - E(I) / RT )`. This weight is added to the 
`overallZ` partition function. Furthermore, we add the weight to a respective
spot associated partition function `spotZ`, if the interaction `I` spans the spot, ie.
the spot's indices are within the interaction subsequences of `I`. If none of
the spots if spanned by `I`, the `noSpotZ` partition function is increased by
`bw(I)`. The final probability of a spot is than given by `spotZ/overallZ` and
the probability of interactions not covering any of the tracked spots is 
computed by `noSpotZ/overallZ` and reported for the pseudo-spot encoding `0&0`
(since indexing starts with 1).

*NOTE* and be aware that the probabilities are *only estimates* since IntaRNA
is not considering (in default prediction mode) all possible interactions due 
to its heuristic (see [discussion about suboptimal interactions](#subopts)).
Nevertheless, since the Boltzmann probabilities are dominated by the low(est)
energy interactions, we consider the probability estimates as meaningful!


<br />
<a name="accessibility" />

### Accessibility and unpaired probabilities

Accessibility describes the availability of an RNA subsequence for intermolecular
base pairing. It can be expressed in terms of the probability of the subsequence
to be unpaired (its *unpaired probability* *Pu*).

A limited accessibility, i.e. a low unpaired probability, can be incorporated into
the RNA-RNA interaction prediction by adding according energy penalties. 
These so called *ED* values are transformed unpaired probabilities, i.e. the
penalty for a subsequence partaking in an interaction is given by *ED=-RT log(Pu)*, 
where *Pu* denotes the unpaired probability of the subsequence. Within the 
IntaRNA energy model, *ED* values for both interacting subsequences are considered.

Accessibility incorporation can be disabled for query or target sequences using
`--qAcc=N` or `--tAcc=N`, respectively.

A setup of `--qAcc=C` or `--tAcc=C` (default) enables accessibility computation 
using the Vienna RNA package routines for query or target sequences, respectively.


<a name="accLocalGlobal" />

#### Local versus global unpaired probabilities

Exact computation of unpaired probabilities (*Pu* terms) is considers all possible
structures the sequence can adopt (the whole structure ensemble). This is referred
to as *global unpaired probabilities* as computed e.g. by **RNAup**.

Since global probability computation is (a) computationally demanding and (b) not
reasonable for long sequences, local RNA folding was suggested, which also enables
according *local unpaired probability* computation, as e.g. done by **RNAplfold**.
Here, a folding window of a defined length 'screens' along the RNA and computes
unpaired probabilities within the window (while only intramolecular base pairs 
within the window are considered).

IntaRNA enables both global as well as local unpaired probability computation.
To this end, the sliding window length has to be specified in order to enable/disable
local folding.

##### Use case examples global/local unpaired probability computation

The use of global or local accessibilities can be defined independently 
for query and target sequences using `--qAccW|L` and `--tAccW|L`, respectively.
Here, `--?AccW` defines the sliding window length (0 sets it to the whole sequence length)
and `--?AccL` defines the maximal length of considered intramolecular base pairs,
i.e. the maximal number of positions enclosed by a base pair
(0 sets it to the whole sequence length). Both can be defined
independently while respecting `AccL <= AccW`.
```bash
# using global accessibilities for query and target
IntaRNA [..] --qAccW=0 --qAccL=0 --tAccW=0 --qAccL=0
# using local accessibilities for target and global for query
IntaRNA [..] --qAccW=0 --qAccL=0 --tAccW=150 --qAccL=100
```


<a name="accConstraints" />

#### Constraints for accessibility computation

For some RNAs additional accessibility information is available. For instance,
it might be known from experiments that some subsequence is unpaired or already
bound by some other factor. The first case (unpaired) makes such regions 
especially interesting for interaction prediction and should result in no ED
penalties for these regions. In the second case (blocked) the region should be 
excluded from interaction prediction.

To incorporate such information, IntaRNA provides the possibility to constrain
the accessibility computation using the `--qAccConstr` and `--tAccConstr` 
parameters. Both take a string encoding for each sequence position whether it is

- `.` unconstrained
- `x` for sure accessible (unpaired)
- `p` paired intramolecularly with some other position of this RNA
- `b` blocked by some other interaction (implies single-strandedness)

Note, *blocked* regions are currently assumed to be bound single-stranded by some
other factor and thus are *treated as unpaired* for ED computation.

```bash
# constraining some central query positions to be blocked by some other molecules
IntaRNA [..] --query="GGGGGGGCCCCCCC" \
        --qAccConstr="...bbbb......."
```

It is also possible to provide a more compact index-range-based encoding of the
constraints, which is especially useful for longer sequences or if you have only
a few constrained regions. To this end, one can provide a comma-separated list 
of index ranges that are prefixed with the according constraint letter from 
above and a colon. Best check the following examples, which should give a good
idea how to use. Note, indexing is supposed to be based on a minimal index of 1
and all positions not covered by the encoding are assumed to be unconstrained
(which must not to be encoded explicitely).

```bash
# applying the same constraints by different encodings to query and target
# example 1
IntaRNA [..] --qAccConstr="...bbbb....." --tAccConstr="b:4-7"
# example 2
IntaRNA [..] --qAccConstr="..bb..xxp.bb" --tAccConstr="b:3-4,11-12,x:7-8,p:9-9"
```



<a name="accFromFile" />

#### Read/write accessibility from/to file or stream

It is possible to read precomputed accessibility values from file or stream to
avoid their runtime demanding computation. To this end, we support the following
formats

| Input format | produced by |
| ---- | --- |
| RNAplfold unpaired probabilities | `RNAplfold -u` or `IntaRNA --out=*Pu:` |
| RNAplfold-styled ED values | `IntaRNA --out=*Acc:` |

The **RNAplfold** format is a table encoding of a banded upper triangular matrix 
with band width l. First row contains a header comment on the data starting with
`#`. Second line encodes the column headers, i.e. the window width per column.
Every successive line starts with the index (starting from 1) of the window end
followed by a tabulator separated list for each windows value in increasing
window length order. That is, column 2 holds values for window length 1, column 
3 for length 2, ... . The following provides a short output/input 
example for a sequence of length 5 with a maximal window length of 3.

```
#unpaired probabilities
 #i$	l=1	2	3	
1	0.9949492	NA	NA	
2	0.9949079	0.9941056	NA	
3	0.9554214	0.9518663	0.9511048		
4	0.9165814	0.9122866	0.9090283		
5	0.998999	0.915609	0.9117766		
6	0.8549929	0.8541667	0.8448852		

```


##### Use case examples for read/write accessibilities and unpaired probabilities

If you have precomputed data, e.g. the file `plfold_lunp` with unpaired probabilities
computed by **RNAplfold**, you can run

```bash
# fill accessibilities from RNAplfold unpaired probabilities
IntaRNA [..] --qAcc=P --qAccFile=plfold_lunp
# fill accessibilities from RNAplfold unpaired probabilities via pipe
cat plfold_lunp | IntaRNA [..] --qAcc=P --qAccFile=STDIN
```

Another option is to store the accessibility data computed by IntaRNA for 
successive calls using 

```bash
# storing and reusing (target) accessibility (Pu) data for successive IntaRNA calls
IntaRNA [..] --out=tPu:intarna.target.pu
IntaRNA [..] --tAcc=P --tAccFile=intarna.target.pu
# piping (target) accessibilities (ED values) between IntaRNA calls
IntaRNA [..] --out=tAcc:STDOUT | IntaRNA [..] --tAcc=E --tAccFile=STDIN
```


Note, for *multiple sequences* in FASTA input, one can also load the
accessibilities (for *all* sequencces) from file. To this end, the file names
have to be prefixed with with `s` and the according sequence's number (where indexing
starts with 1) within the FASTA input using a common suffix after the index.
This suffix is to be provided to the according `--?AccFile` argument. 
The files generated by `--out=?Acc:...` are already conform to this requirement,
such that you can use the use case examples from above also for multi-sequence
FASTA input. 
Note, this is not supported for a piped setup (e.g. via `--out=tAcc:STDOUT`
as shown above), since this does not produce the according output files!




<br /><br />
<a name="multithreading" />

## Multi-threading and parallelized computation

IntaRNA supports the parallelization of the target-query-combination processing. 
The maximal number of threads to be used can be specified using the `--threads` parameter.
If `--threads=k != 1`, than *k* predictions are processed in parallel. A value of
0 requests the maximally available number of threads for this machine.

When using parallelization, you should have the following in mind:

- The memory consumption will be (much) higher, since each thread runs an independent
  prediction (with according memory consumption). Thus, ensure you have enough
  RAM available when using many threads of memory-demanding 
  [prediction modes](#predModes).
 
The support for multi-threading can be completely disabled before compilation
using `configure --disable-multithreading`.







<br /><br /><br /><br />
<a name="lib" />

# Library for integration in external tools

The IntaRNA package also comes with a C++ library `libIntaRNA.a` containing the core classes
and functionalities used within the IntaRNA tool. The whole library comes with
an `IntaRNA` namespace and exhaustive class and member API documentation that is
processed using doxygen to generate html/pdf versions.

When IntaRNA is build while `pkg-config` is present, according pkg-config
information is generated and installed too.

## Mandatory `Easylogging++` initalization !

Since IntaRNA makes heavy use of the `Easylogging++` library, you have to add (and adapt) 
the following code to your central code that includes the `main()` function:
```[c++]
// get central IntaRNA-lib definitions and includes
#include <IntaRNA/general.h>
// initialize logging for binary
INITIALIZE_EASYLOGGINGPP

[...]

int main(int argc, char **argv){

[...]
		// set overall logging style
		el::Loggers::reconfigureAllLoggers(el::ConfigurationType::Format, std::string("# %level : %msg"));
		// no log file output
		el::Loggers::reconfigureAllLoggers(el::ConfigurationType::ToFile, std::string("false"));
		el::Loggers::reconfigureAllLoggers(el::ConfigurationType::ToStandardOutput, std::string("true"));
		// set additional logging flags
		el::Loggers::addFlag(el::LoggingFlag::DisableApplicationAbortOnFatalLog);
		el::Loggers::addFlag(el::LoggingFlag::LogDetailedCrashReason);
		el::Loggers::addFlag(el::LoggingFlag::AllowVerboseIfModuleNotSpecified);

		// setup logging with given parameters
		START_EASYLOGGINGPP(argc, argv);
[...]
}
```

Note further, to get the library correctly working the following compiler 
flags are used within the IntaRNA configuration:
```[bash]
    CXXFLAGS=" -DELPP_FEATURE_PERFORMANCE_TRACKING -DELPP_NO_DEFAULT_LOG_FILE "
``` 



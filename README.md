[![Build Status](https://travis-ci.org/BackofenLab/IntaRNA.svg?branch=master)](https://travis-ci.org/BackofenLab/IntaRNA)

# IntaRNA version 2.*

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
- [IntaRNA: efficient prediction of bacterial sRNA targets incorporating target site accessibility and seed regions](http://dx.doi.org/10.1093/bioinformatics/btn544)
  Anke Busch, Andreas S. Richter, and Rolf Backofen, 
  Bioinformatics, 24 no. 24 pp. 2849-56, 2008, DOI(10.1093/bioinformatics/btn544).
- [CopraRNA and IntaRNA: predicting small RNA targets, networks and interaction domains](http://dx.doi.org/10.1093/nar/gku359)
  Patrick R. Wright, Jens Georg, Martin Mann, Dragos A. Sorescu, Andreas S. Richter, Steffen Lott, Robert Kleinkauf, Wolfgang R. Hess, and Rolf Backofen
  Nucleic Acids Research, 42 (W1), W119-W123, 2014, DOI(10.1093/nar/gku359).



<br /><br /><br /><br />
<a name="doc" />
# Documentation

## Overview

The following topics are covered by this documentation:

- [Installation](#install)
  - [IntaRNA via conda](#instconda)
  - [Dependencies](#deps)
  - [Cloning from github](#instgithub)
  - [Source code distribution](#instsource)
  - [IntaRNA docker container](#instdocker)
- [Usage and Parameters](#usage)
  - [Just run ...](#defaultRun)
  - [Prediction modes, their features and emulated tools](#predModes)
  - [Interaction restrictions](#interConstr)
  - [Seed constraints](#seed)
  - [Output modes](#outmodes)
  - [Suboptimal RNA-RNA interaction prediction and output restrictions](#subopts)
  - [Energy parameters and temperature](#energy)
  - [Additional output files](#outFiles)
    - [Minimal energy profiles](#profileMinE)
    - [Accessibility and unpaired probabilities](#accessibility)
      - [Local versus global unpaired probabilities](#accLocalGlobal)
      - [Read/write accessibility from/to file or stream](#accFromFile)
  - [Multi-threading and parallelized computation](#multithreading)



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
conda install intarna
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
  (ensure the following libraries are installed; or install all e.g. in Ubuntu via package `libboost-all-dev`)
  - libboost_regex
  - libboost_program_options
  - libboost_filesystem
  - libboost_system
- [Vienna RNA package](http://www.tbi.univie.ac.at/RNA/) version >= 2.3.0
- if [cloning from github](#instgithub): GNU autotools (automake, autoconf, ..)

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
- `--with-RNA` : the prefix where the Vienna RNA package is installed
- `--with-boost` : the prefix where the boost library is installed


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




<br /><br /><br /><br />
<a name="usage" />
# Usage and parameters

IntaRNA comes with a vast variety of ways to tune or enhance *YOUR* RNA-RNA prediction.
To this end, different [prediction modes](#predModes) are implemented that allow
to balance predication quality and runtime requirement. Furthermore, it is 
possible to define 
[interaction restrictions](#interConstr),
[seed constraints](#seed), 
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
minimal interaction energy without consideratio whether or not the interacting 
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

Furthermore, the region where interactions are supposed to occur can be restricted
for target and query independently. To this end, a list of according index pairs
can be provided using `--qRegion` and `--tRegion`, respectively. The indexing 
starts with 1 and should be in the format `from1-end1,from2-end2,..` using
integers.

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

Seed constraint usage can be globally disabled using the `--noSeed` flag.






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
interaction output in CSV format (using `;` as separator).

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

- `id1` : id of first sequence
- `id2` : id of second sequence
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




<br /><br />
<a name="subopts" />
## Suboptimal RNA-RNA interaction prediction and output restrictions

Besides the identification of the optimal (e.g. minimum-free-energy) RNA-RNA 
interaction, IntaRNA enables the enumeration of suboptimal interactions. To this
end, the argument `-n N` or `--outNumber=N` can be used to generate up to `N`
interactions for each query-target pair (including the optimal one). Note, the
suboptimal enumeration is increasingly sorted by energy.

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
- `qAcc:`|`tAcc:` the [query/target's ED accessibility values](#accessibility) (RNAplfold-like format), respectively
- `qPu:`|`tPu:` the [query/target's unpaired probabilities](#accessibility) (RNAplfold format), respectively




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
plot( d[,1], d[,2], xlab="sequence index", ylab="minimal energy", type="l", col="blue", lwd=2)
abline(h=0, col="red", lty=2, lwd=2)
```

![Minimal interaction energy profile of an RNA](/doc/figures/profile-minE.png?raw=true "Minimal interaction energy profile of an RNA")

This plot reveals two less but still stable (*E* below 0) interaction sites beside the
mfe interaction close to the 5'-end of the molecule.


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
image(1:nrow(d),1:ncol(t),t, col = heat.colors(100), xlab="sequence 1 index", ylab="sequence 2 index");
box();
```

The following plot (for the [minimal energy profile](#profileMinE) example from
above) reveals, that the alternative stable (*E*<0) interactions all involve the
mfe-site in the second sequence and are thus less likely to occure. 

![Minimal interaction energy index pair information](/doc/figures/pair-minE.png?raw=true "Minimal interaction energy index pair information")




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






<br /><br />
<a name="multithreading" />
## Multi-threading and parallelized computation

IntaRNA supports the parallelization of the target-query-combination processing. 
The maximal number of threads to be used can be specified using the `--threads` parameter.
If `--threads=k > 0`, than *k* predictions are processed in parallel.

When using parallelization, you should have the following things in mind:

- Most of the IntaRNA runtime (in heuristic prediction mode) 
  is consumed by [accessibility computation](#accessibility) 
  (if not [loaded from file](#accFromFile)). 
  Currently, due to some thread-safety issues with the 
  routines from the Vienna RNA package, the IntaRNA
  accessibility computation is done serially. This significantly reduces the
  multi-threading effect when running IntaRNA in the fast heuristic mode (`--mode=H`).
  If you run a non-heuristic prediction mode, multi-threading will show a more
  dramatic decrease in runtime performance, since here the interaction prediction
  is the computationally more demanding step.
- The memory consumption will be much higher, since each thread runs an independent
  prediction (with according memory consumption). Thus, ensure you have enough
  RAM available when using many threads of memory-demanding 
  [prediction modes](#predModes).
 
The support for multi-threading can be completely disabled before compilation
using `configure --disable-multithreading`.

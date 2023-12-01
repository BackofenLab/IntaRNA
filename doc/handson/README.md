
# Hands-on examples for the application of IntaRNA

In the following, we list some examples how to use IntaRNA in specific application scenarios.
The examples are presented in detail in our publication

- [How to do RNA-RNA interaction prediction? A use-case driven
handbook using IntaRNA](http://www.bioinf.uni-freiburg.de/Subpages/publications.html?de#Raden-IntaRNA-handson.abstract)
  - Martin Raden and Milad Miladi
  - Springer (in press, DOI to come)




## Example 2.1: Local Installation via `conda` (Linux System with `bash` Shell)

First, we need to install `conda` (if not already available).
To that end, we assume you use a Linux-like system with a `bash` command line shell.
On Microsoft Windows, you can emply the Windows Subsystem for Linux (WSL) for that purpose.

```sh
# download Miniconda installer
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
# note: ’-b’ runs the installer in non-interactive silent mode (agreeing to licence etc.)
# starting conda installation process
bash ~/miniconda.sh -b -p $HOME/miniconda
```

Everytime when starting a new shell (and thus before running IntaRNA), we need to ensure `conda` is among the tools we can run, i.e. its tools are among the program search path.
You can simplify this, by adding the `export` line to your `.bashrc` script in your home directory.

```sh
# add conda to your program search path
export PATH="$HOME/miniconda/bin:$PATH"
# installing IntaRNA via conda
conda install -c conda-forge -c bioconda "intarna>3" python=3.10 # python needed for viennarna package :/
# checking installed IntaRNA version
IntaRNA --version
```

It might be that your `conda` environment comes with a preinstalled `python` version that exceeds the dependency specifications of the `viennarna` package, a central dependency of IntaRNA.
In that case, you need to specify a `python` version that meets the requirements within the `conda install` call.
For example at the time of writing (Dec 2023), I used for `conda` 23.10.0 to install IntaRNA 3.3.2 with `viennarna` 2.5.1 the following call.

```sh
conda install -c conda-forge -c bioconda "intarna>3" python=3.10 # needed for viennarna 2.5.1
```


## Example 3.1: Standard Input and Output

The bacterial sRNA OxyS of Echerichia coli is known to interact with the mRNA fhlA as shown by [Argaman and Altuvia (2000)](https://doi.org/10.1006/jmbi.2000.3942).
From that publication, we extracted the sequences shown in Figure 7 provided here in the files [fhlA.fasta](fhlA.fasta) and [OxyS.fasta](OxyS.fasta).
We use `--tIdxPos0` to index the mRNA in relation to the position of its start codon.    

```sh
# minimal call with with explicit fhlA index information
IntaRNA -t fhlA.fasta -q OxyS.fasta --tIdxPos0=-53
```

produces 

```
fhlA|NC000913|-53..+60|doi:10.1006/jmbi.2000.3942:Fig7
           -15                      +7
             |                      |
5'-AGU...GCUU       A-- AA     UG    AUAC...GAC-3'
             UCCUGGA   C  CAAAA  UCAU
             +++++++   |  |||||  |||:
             AGGACCU   G  GUUUU  AGUG
      3'-GCCU       CUA GC     CA    CAAC...AAG-5'
             |                      |
           104                      81
OxyS|NC_000913|56..164|doi:10.1006/jmbi.2000.3942:Fig7

interaction energy = -5.59 kcal/mol
```

Therein, base pairs denoted with `+` are part of putative seed regions (typically most reliable) while `|` and `:` represent `AU,CG` and `GU` base pairs, respectively.



## Example 3.2: Detailed Output of Exact Prediction

To investigate the details of interactions, we suggest the `IntaRNAexact` personality.

```sh
IntaRNAexact --intLenMax=60 -t fhlA.fasta --tIdxPos0=-53 -q OxyS.fasta --outMode=D
```

`--outMode=D` provides more details for the interaction:

```

fhlA|NC000913|-53..+60|doi:10.1006/jmbi.2000.3942:Fig7
           -15     -9
             |     |
5'-AGU...GCUU       ACAA...GAC-3'
             UCCUGGA
             +++++++
             AGGACCU
      3'-GCCU       CUAG...AAG-5'
             |     |
           104     98
OxyS|NC_000913|56..164|doi:10.1006/jmbi.2000.3942:Fig7

interaction seq1   = -15..-9
interaction seq2   = 98..104

interaction energy = -5.57 kcal/mol
  = E(init)        = 4.1
  + E(loops)       = -15.6
  + E(dangleLeft)  = -0.4
  + E(dangleRight) = -0.96
  + E(endLeft)     = 0.5
  + E(endRight)    = 0.5
    : E(hybrid)    = -11.86
  + ED(seq1)       = 0.62
    : Pu(seq1)     = 0.36569
  + ED(seq2)       = 5.67
    : Pu(seq2)     = 0.000101064

seed seq1   = -15..-9
seed seq2   = 98..104
seed energy = -5.57
seed ED1    = 0.62
seed ED2    = 5.67
seed Pu1    = 0.36569
seed Pu2    = 0.000101064
```

It is also useful to study suboptimal interactions using `-n`:

```sh
IntaRNAexact --intLenMax=60 -t fhlA.fasta --tIdxPos0=-53 -q OxyS.fasta \
  -n 10 --outMode=C --outCsvSort=E --outCsvCols="E,Pu2,subseqDB,hybridDB"
```

Note: `\` at the line end will continue the command in the next line.
You can remove the `\` along with the line break if needed.

```
E;Pu2;subseqDB;hybridDB
-5.57;0.000101064;-15UCCUGGA&98UCCAGGA;-15|||||||&98|||||||
-5.54;5.61948e-07;-15UCCUGGAACAACAAAAUGUCAU&81GUGAACUUUUGCGGAUCUCCAGGA;-15|||||||.|..|||||..||||&81||||..|||||..|...|||||||
-3.98;5.61948e-07;-15UCCUGGAACAACAAAAUGUCA&82UGAACUUUUGCGGAUCUCCAGGA;-15|||||||.|..|||||..|||&82|||..|||||..|...|||||||
-3.87;5.52904e-07;-15UCCUGGAACAACAAAAUGUCAUAU&79ACGUGAACUUUUGCGGAUCUCCAGGA;-15|||||||.|..|||||..||||.|&79|.||||..|||||..|...|||||||
-3.62;0.0167597;34CAAGGGU&24ACCCUUG;34|||||||&24|||||||
-3.42;6.6094e-07;-15UCCUGGAACAACAAAA&87UUUUGCGGAUCUCCAGGA;-15|||||||.|..|||||&87|||||..|...|||||||
-3.22;5.99629e-07;-15UCCUGGAACAACAAAAUGUC&83GAACUUUUGCGGAUCUCCAGGA;-15|||||||.|..|||||..||&83||..|||||..|...|||||||
-2.74;5.44005e-07;-15UCCUGGAACAACAAAAUGUCAUAUACAC&75GCCAACGUGAACUUUUGCGGAUCUCCAGGA;-15|||||||.|..|||||..||||.....|&75|.....||||..|||||..|...|||||||
-2.72;6.71751e-07;-15UCCUGGAACAACAAA&88UUUGCGGAUCUCCAGGA;-15|||||||.|..||||&88||||..|...|||||||
-2.7;6.6094e-07;-15UCCUGGAACAACAAAAUGU&85ACUUUUGCGGAUCUCCAGGA;-15|||||||.|..|||||.||&85|||||||..|...|||||||
```

The output reveals the second interaction site on rank 5 that was reported by [Argaman and Altuvia (2000)](https://doi.org/10.1006/jmbi.2000.3942).


## Example 3.3:  Corrected Prediction of 2nd RRI in Two-site Concurrent Model

Given the first interaction site predicted in Example 3.2, we can predict an interaction model of two concurrently formed interaction sites using the following command.

```sh
IntaRNAexact --intLenMax=60 -t fhlA.fasta --tIdxPos0=-53 -q OxyS.fasta \
  --tAccConstr="b:-15--9" --qAccConstr="b:98-104" --energyAdd="-5.57"
```

`--energyAdd` provides the interaction energy of the first site while the `--*AccConstr` accessibility constraints ensure that the regions involved in the first interaction are considered blocked for the prediction of the second interaction site.
Thus, IntaRNA predicts a concurrent second site (in accordance with literature) and reports the energy of the join interaction formation of both sites.

```
fhlA|NC000913|-53..+60|doi:10.1006/jmbi.2000.3942:Fig7
           +34     +40
             |     |
5'-AGU...ACAA       UGUU...GAC-3'
             CAAGGGU
             +++++++
             GUUCCCA
3'-GCC...UGAA       AUUU...AAG-5'
             |     |
            30     24
OxyS|NC_000913|56..164|doi:10.1006/jmbi.2000.3942:Fig7

interaction energy = -9.19 kcal/mol
```


## Example 3.4: Anchored Prediction

When specific base pairs of an interaction are known, e.g. from mutation experiments, this information can be provided to constrain the IntaRNA prediction.
Here we use the interaction of the bacterial sRNA GcvB with its mRNA target phoB.
Sequences [GcvB.fasta](GcvB.fasta) and [phoB.fasta](phoB.fasta) were taken from the NCBI reference genome NC_000913.

```sh
IntaRNAexact -q GcvB.fasta -t phoB.fasta --tIdxPos0=-200 --intLenMax=60 -n 2
```

When using an informed standard prediction mode, the following interactions are predicted.

```

phoB|NC_000913|b1130|-200..+100|genom-subsequence
           -56                     -34
             |                     |
5'-GAG...CCCC      C     AUC        CUAU...GCC-3'
             CAUAAC ACAUA   GCGUUACA
             ||:||| |||:|   ++++++++
             GUGUUG UGUGU   UGUAGUGU
3'-UUU...GUUU      U     ---        UGGC...UCA-5'
             |                     |
            85                     66
GcvB|NC_000913

interaction energy = -9.44 kcal/mol

phoB|NC_000913|b1130|-200..+100|genom-subsequence
           -18      -11
             |      |
5'-GAG...AUUA        AGAA...GCC-3'
             AGACAGGG
             ++++++++
             UCUGUCCC
3'-UUU...CCUG        AUUU...UCA-5'
             |      |
           159      152
GcvB|NC_000913

interaction energy = -9.34 kcal/mol
```

Here the known interaction from [(Coornaer et al., 2013)](https://doi.org/10.1371/journal.pgen.1003156) is only second in rank.
Incorporating the experimental information from the publication using

```sh
IntaRNAexact -q GcvB.fasta -t phoB.fasta --tIdxPos0=-200 --intLenMax=60 -n 2 \
  --seedTQ="-17|||||&154|||||"
```

provides variants of the validated interaction site.

```
phoB|NC_000913|b1130|-200..+100|genom-subsequence
           -18      -11
             |      |
5'-GAG...AUUA        AGAA...GCC-3'
             AGACAGGG
             |+++++||
             UCUGUCCC
3'-UUU...CCUG        AUUU...UCA-5'
             |      |
           159      152
GcvB|NC_000913

interaction energy = -9.34 kcal/mol

phoB|NC_000913|b1130|-200..+100|genom-subsequence
           -17     -11
             |     |
5'-GAG...UUAA       AGAA...GCC-3'
             GACAGGG
             +++++||
             CUGUCCC
3'-UUU...CUGU       AUUU...UCA-5'
             |     |
           158     152
GcvB|NC_000913

interaction energy = -8.53 kcal/mol
```



## Example 3.5: Seed-region-constraint Prediction

In case seed regions of RNAs, i.e. regions known to be involved in interactions, are known, their specification can provide a powerful tool to enhance the prediction quality.
For instance, GcvB is known to feature three regions R1, R2 and R3 that can be involved in regulatory interactions.

Regions R1 and R2 correspond to the positions ranges 75-93 and 129-145 in the GcvB sequence of S. Typhimurium provided in [GcvB.ST.fasta](GcvB.ST.fasta).
When predicting interactions with its the regulatory region of its mRNA target ilvE provided in [ilvE.fasta](ilvE.fasta) using

```sh
IntaRNAexact -q GcvB.ST.fasta -t ilvE.fasta --tIdxPos0=-200 --intLenMax=60
``` 

the top-ranked interaction is heavily disrupting local intra-molecular structure elements.

```
ilvE|NC_003197|STM3903|-200..+100|genom-subsequence
           -56                          -32
             |                          |
5'-GGU...GAUC   ---  G    C   C          AAAU...GGU-3'
             UGC   CA AGCG UGC ACAUCACAAC
             |||   || |:|: ::| ++++++++++
             ACG   GU UUGU GUG UGUAGUGUUG
3'-UUC...GUUA   UUU  G    U   U          GCAU...UCA-5'
             |                          |
            91                          64
GcvB|NC_003197

interaction energy = -8.32 kcal/mol
```


When restricting seed prediction to the known regions R1 and R2, using

```sh
IntaRNAexact -q GcvB.ST.fasta -t ilvE.fasta --tIdxPos0=-200 --intLenMax=60 \
  --seedQRange="75-93,129-145"
```

the correct and experimentally verified interaction is predicted.

```
ilvE|NC_003197|STM3903|-200..+100|genom-subsequence
           -37      -30
             |      |
5'-GGU...ACAU        AUCC...GGU-3'
             CACAACAA
             ++++++++
             GUGUUGUU
3'-UUC...GUUU        GUGU...UCA-5'
             |      |
            85      78
GcvB|NC_003197

interaction energy = -6.38 kcal/mol
```


## Example 3.6: Target Screen on Genomic Regions

Using the regulatory regions around all start codons of mRNAs in Echerichia coli provided in [NC_000913.fa](NC_000913.fa), we can do target prediction for the sRNA ChiX given in [ChiX_NC_000913.fasta](ChiX_NC_000913.fasta).

To prepare multiple target screens (not done here, but typically the case), we first generate and store the accessibility information of all putative targets.

```sh
# Precomputation of target accessibility via dummy call (pseudo query + no output)
IntaRNA -t NC_000913.fa --out=tAcc:NC_000913.ed -q AAA -n 0 --out=/dev/null --noSeed
```

Afterwards, we can do a fast target scan using

```sh
# Multi-threaded target screen using precomputed target accessibility data
IntaRNA -t NC_000913.fa -q ChiX_NC_000913.fasta \
 --intLenMax=60 \
 --outMode=C --outCsvCols="E,id1,start1,end1,id2,start2,end2" --outCsvSort=E \
 --tAcc=E --tAccFile=NC_000913.ed \
 --threads=0 \
 > NC_000913.ChiX.csv
```

and the resulting top-5 predictions are

```
E;id1;start1;end1;id2;start2;end2
-66.23;b0481;20;78;ChiX_NC_000913;2;60
-16.41;b0681;182;193;ChiX_NC_000913;45;56
-15.45;b0619;164;175;ChiX_NC_000913;46;57
-12.9;b0352;157;192;ChiX_NC_000913;45;84
-12.51;b0847;224;253;ChiX_NC_000913;26;55
```

Note, `b0481` is the gene encoding ChiX and thus highly complementary.
This shows that intra-molecular structure prediction and interaction prediction are quite related tasks.
  
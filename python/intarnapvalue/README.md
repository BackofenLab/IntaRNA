# IntaRNApvalue
A tool for calculating p-values of IntaRNA energy scores.

### How it works:
This tool shuffles one or more of the given sequences (depending on shuffle mode) while preserving mono- and dinucleotide frequences.
Thus the resulting sequences are very similar with similar properties to the original one.
It then calculates the IntaRNA interaction scores for these newly generated sequences and uses them to approximate the p-value of the original target/query combination.
This p-value is a score for the likelihood that a sequence with a better interaction score than the original ones can randomly occur.
It can thus be used as a measurement for how good an interaction between two sequences is.

## Dependencies:
- IntaRNA
##### And if you want to run it with python or compile it yourself:
- Python 3 (>= 3.6)
- numpy
- scipy

## Installation:
You can either run this tool directly with python, install it with setuptools or run it as a compiled binary.
##### Run with python directly:
You don't have to install it, if you just plan on running it from a certain directory.
Simply copy the "intarnapvalue" folder where you want to run it from.
You have to get the dependencies yourself in this case.

To use the tool from anywhere however, you need to install it as a python module.
All dependencies (except IntaRNA) will be installed automatically.
This way it can also be installed in a virtual environment.
```console
python setup.py install
```

##### Run as binary:
Go into the bin directory. You can find pre-compiled builds for linux and windows for x64 directly.
You don't need to install any dependencies for running these, on Windows however you will need Microsoft Visual C++ Redistributable.
If you want to compile it yourself, you need to get the dependencies and PyInstaller.
I had problems with numpy 1.17.0, so I had to use 1.16.4.
Simply run the build.py with the python binary that has the dependencies and PyInstaller installed.
You will find your binary in the build folder.

## Usage:
Run it as a module without installing (you need to be in the parent dir of intarnapvalue):
```console
python3 -m intarnapvalue <arguments>
```

Run it as a module from anywhere with python (need to run setup.py first, see [Installation](#installation)):
```console
python3 -m intarnapvalue <arguments>
```
Or use as a compiled binary:
```console
./IntaRNApvalue <arguments>
```

Or import it and use it from within other python code:
```python
from intarnapvalue.intarna_pvalue import IntaRNApvalue
IntaRNApvalue(['--flag1', 'arg1'])
```

## Arguments:

| Flag                | Value                  | Default | Description          |
| ------------------  |:---------------------- | :------ | -------------------- |
| -h, --help          |                        |         | Gives detailed command help.  |
| -q, --query         | Sequence or FASTA file |         | Query as a raw sequence or a file in FASTA format. Takes the first sequence from the file. |
| -t, --target        | Sequence or FASTA file |         | Target as a raw sequence or a file in FASTA format. Takes the first sequence from the file. |
| -c, --cardinality   | Integer                |         | How many sequence pairs are randomly permuted and considered for p-value calculation. |
| -m, --shuffle-mode  | {q, t, b}              |         | Which sequence will be shuffled: Query, Target or both. |
| -d, --distribution  | {gev, gumbel, gauss}   | gev     | (optional) The distribution used for p-value calculation: Generalized Extreme Value Distribution, Gumbel Distribution or Gauss. |
| -o, --output        | {pvalue, scores}       | pvalue  | (optional) If set to p-value, will only output p-value. If set to scores, will output every score from randomly generated sequences, but no p-value. |
| --threads           | 0 - {max threads}      | 0       | (optional) How many threads IntaRNA uses for score calculation. If set to 0 it will use all available threads. |
| --randSeed          | Integer                | None    | (optional) The seed used for generating random sequences. |
| -p, --parameterFile | file name              | None    | (optional) parameter file to be used in IntaRNA calls to further guide predictions. |

## Example
An example use-case could look like this with sequences in raw format:
```console
$ python3 -m intarnapvalue --query GCUGAAAAACAUAACCCAUAAAAUGCUAGCUGUACCAGGAACCA --target GGUUUCUUCGCCUCUGCGUUCACCAAAGUGUUCACCC --scores 10000 --shuffle-mode b --threads 0
0.2421277140114205
```
Or you could specify any of the sequences as files:
```console
$ python3 -m intarnapvalue --query query.fasta --target target.fasta --scores 10000 --shuffle-mode b --threads 0
0.2421277140114205
```
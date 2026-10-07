# IntaRNA C++ API

IntaRNA predicts interactions between two RNA molecules while accounting for
the energy needed to make their binding sites accessible and for configurable
seed interactions. The command-line tool and the `libIntaRNA` C++ library share
the same prediction machinery.

This reference describes the library in the `IntaRNA` namespace. For installation,
command-line options, examples, and citations, see the
[IntaRNA user guide](https://backofenlab.github.io/IntaRNA/).
The online API follows the development branch (`master`); generate the
documentation from a release checkout when working with that release.

## Finding your way around

| Task | Starting points |
| --- | --- |
| Represent sequences and interaction sites | @ref IntaRNA::RnaSequence, @ref IntaRNA::IndexRange, @ref IntaRNA::Interaction |
| Compute or load accessibility penalties | @ref IntaRNA::Accessibility, @ref IntaRNA::AccessibilityVrna, @ref IntaRNA::AccessibilityFromStream |
| Evaluate interaction energies | @ref IntaRNA::InteractionEnergy, @ref IntaRNA::InteractionEnergyVrna |
| Predict minimum-energy interactions | @ref IntaRNA::Predictor, @ref IntaRNA::PredictorMfe, @ref IntaRNA::PredictorMfe2dHeuristic |
| Work with ensemble-based predictions | @ref IntaRNA::PredictorMfeEns, @ref IntaRNA::PredictorMfeEns2d |
| Specify seeds and helices | @ref IntaRNA::SeedConstraint, @ref IntaRNA::SeedHandler, @ref IntaRNA::HelixConstraint, @ref IntaRNA::HelixHandler |
| Collect predictions or track computations | @ref IntaRNA::OutputHandler, @ref IntaRNA::OutputConstraint, @ref IntaRNA::PredictionTracker |

Use the **Classes** menu for the class list and inheritance hierarchy, **Files**
for public headers, or the search box for a particular symbol. Concrete classes
document their supported constraints and algorithmic tradeoffs.

## How the pieces fit together

1. Create an @ref IntaRNA::RnaSequence for each RNA and choose an
   @ref IntaRNA::Accessibility implementation for each sequence. Accessibility
   describes the energy penalty for making a subsequence available for binding.
2. Wrap the second sequence's accessibility in @ref IntaRNA::ReverseAccessibility
   and pass both accessibility objects to an @ref IntaRNA::InteractionEnergy
   implementation. The reversed view lets the energy model handle antiparallel
   pairing consistently.
3. Select a concrete @ref IntaRNA::Predictor and provide the energy model, an
   @ref IntaRNA::OutputHandler, and any required seed or helix handlers. Calling
   @ref IntaRNA::Predictor::predict delivers predicted interactions to the output
   handler; optional @ref IntaRNA::PredictionTracker implementations collect
   additional information during prediction.

Many constructors retain references to their inputs. Keep those objects alive
for the lifetime of the objects using them, and consult constructor documentation
for ownership exceptions (including prediction trackers).

## Coordinates and energy units

- Internal sequence positions are zero-based, and @ref IntaRNA::IndexRange uses
  inclusive bounds. Input/output numbering can differ; use the conversions
  supplied by @ref IntaRNA::RnaSequence.
- Sequence 2 is reversed in the energy model. Use
  @ref IntaRNA::ReverseAccessibility and the energy model's conversion methods
  when moving between internal coordinates and reported interactions.
- @ref IntaRNA::E_type stores energies in hundredths of kcal/mol. The conversion
  macros `E_2_Ekcal` and `Ekcal_2_E` in @ref general.h convert between internal
  values and kcal/mol. Preserve the infinity sentinels and use the comparison
  helpers in that header; partition-function values use @ref IntaRNA::Z_type.

## Using the library

Public headers require C++23. The installed `IntaRNA` pkg-config package supplies
include and linker flags; the consuming project must select its C++ language
standard. Follow the
[library integration guide](https://backofenlab.github.io/IntaRNA/#lib), including
the required Easylogging++ initialization, before calling the library.

To build this reference without compiling IntaRNA, run `bash doc/build-api.sh`
from a source checkout with Doxygen and Graphviz installed, then open
`doxygen-doc/html/index.html`. See the
[documentation build guide](https://github.com/BackofenLab/IntaRNA/blob/master/doc/api-documentation.md)
for the Autotools target and publishing setup.

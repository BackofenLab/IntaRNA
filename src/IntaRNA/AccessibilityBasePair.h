
#ifndef INTARNA_ACCESSIBILITYBASEPAIR_H_
#define INTARNA_ACCESSIBILITYBASEPAIR_H_

#include "IntaRNA/general.h"
#include "IntaRNA/RnaSequence.h"
#include "IntaRNA/Accessibility.h"
#include "IntaRNA/AccessibilityConstraint.h"
#include "IntaRNA/NussinovHandler.h"


namespace IntaRNA {

/**
 * Provides accessibility penalties for a Nussinov-like energy scoring where
 * each base pair is scored by #basePairEnergy independently from its structural
 * context.
 *
 * This implementation is mainly for teaching and testing purpose.
 *
 * @author Mostafa Mahmoud
 * @author Martin Mann
 *
 */
class AccessibilityBasePair: public Accessibility {

public:
  /***
   * Constructor of AccessibilityBasePair
   * @param seq The sequence the accessibility data belongs to
   * @param maxLength the maximal length of accessible regions (>0) to be
   *          considered. 0 defaults to the full sequence's length, otherwise
   *          is is internally set to min(maxLength,seq.length).
   * @param accConstr optional accessibility constraint
   * @param basePairEnergy The energy value of the base pairs
   * @param RT The temperature energy constant
   * @param minLoopLength the minimum loop length
   */
  AccessibilityBasePair(
      const RnaSequence& seq,
      const size_t maxLength,
      const AccessibilityConstraint * const accConstr,
      const E_type basePairEnergy = Ekcal_2_E(-1.0),
      const Z_type RT = 1,
      const size_t minLoopLength = 3);

  /***
   * Destructor of AccessibilityBasePair
   */
  virtual ~AccessibilityBasePair();

  /***
   * Returns the accessibility energy value between the indices (from, to)
   * in the given sequence.
   * @param from The start index of the region
   * @param to The end index of the region
   *
   * @return The ED Value
   *
   * @throw std::runtime_error in case it does not hold 0 <= from <= to < seq.length
   */
	virtual E_type getED( const size_t from, const size_t to ) const;

protected:

	//! energy of an individual base pair
	const E_type basePairEnergy;
	//! temperature constant for normalization
	const Z_type RT;
  //! Boltzmann Energy weight
  const Z_type basePairWeight;
  //! minimum length of loops
  const size_t minLoopLength;

  /***
   * Results of getED lookup table
   */
  NussinovHandler::E2dMatrix logPu;

};

}  // namespace IntaRNA

#endif /* INTARNA_ACCESSIBILITYBASEPAIR_H_ */

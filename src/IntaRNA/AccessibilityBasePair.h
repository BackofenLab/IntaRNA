
#ifndef INTARNA_ACCESSIBILITYBASEPAIR_H_
#define INTARNA_ACCESSIBILITYBASEPAIR_H_

#include "general.h"
#include "RnaSequence.h"
#include "Accessibility.h"
#include "AccessibilityConstraint.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


namespace IntaRNA {


class AccessibilityBasePair: public Accessibility {

public:
typedef double P_type;  // Probability type
typedef boost::numeric::ublas::matrix<P_type> P2dMatrix;  // Probability matrix
typedef boost::numeric::ublas::matrix<E_type> E2dMatrix;  // Energy matrix

  /***
   * Constructor of AccessibilityBasePair
	 * @param seq The sequence the accessibility data belongs to
	 * @param maxLength the maximal length of accessible regions (>0) to be
	 *          considered. 0 defaults to the full sequence's length, otherwise
	 *          is is internally set to min(maxLength,seq.length).
	 * @param accConstr optional accessibility constraint
   * @param basePairEnergy The energy value of the base pairs
   * @param RT The temperature energy constant
   */
  AccessibilityBasePair(const RnaSequence& seq, const size_t maxLength,
      const AccessibilityConstraint * const accConstr,
      const E_type basePairEnergy = -1, const E_type RT = 1);

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

private:

  const E_type basePairEnergy;
  const E_type RT;
  const E_type boltzmann = std::exp(-basePairEnergy / RT);
  const size_t minLoopLength = 3;
  const size_t N;

  /***
   * Results of getED lookup table
   */
  E2dMatrix logPu;

  /***
   * Get the partition function Q between the indices (from, to)
   * @param from The start index of the region
   * @param to The end index of the region
   * @param Q Lookup table for Q
   * @param Qb Lookup table for Qb
   * @returns the partition function value
   */
  E_type getQ(const size_t from, const size_t to,
      E2dMatrix &Q, E2dMatrix &Qb);

  /***
   * Get the base partition function Qb between the indices (from, to)
   * @param from The start index of the region
   * @param to The end index of the region
   * @param Q Lookup table for Q
   * @param Qb Lookup table for Qb
   * @returns the partition function value
   */
  E_type getQb(const size_t from, const size_t to,
      E2dMatrix &Q, E2dMatrix &Qb);

  /***
   * Get the base-pair probbility between the indices (from, to)
   * @param from The start index of the region
   * @param to The end index of the region
   * @param Q Lookup table for Q
   * @param Qb Lookup table for Qb
   * @param Pbp Lookup table for base-pair probabilities
   * @returns the probability value
   */
  P_type getPbp(const size_t from, const size_t to,
      E2dMatrix &Q, E2dMatrix &Qb, P2dMatrix &Ppb);

  /***
   * Get the unpaired probbility between the indices (from, to)
   * @param from The start index of the region
   * @param to The end index of the region
   * @param Q Lookup table for Q
   * @param Qb Lookup table for Qb
   * @param Pbp Lookup table for base-pair probabilities
   * @param Pu Lookup table for unpaired probabilities
   * @returns the probability value
   */
  P_type getPu(const size_t from, const size_t to, E2dMatrix &Q, E2dMatrix &Qb,
      P2dMatrix &Ppb, P2dMatrix &Pu);

};

}  // namespace IntaRNA

#endif /* INTARNA_ACCESSIBILITYBASEPAIR_H_ */ 

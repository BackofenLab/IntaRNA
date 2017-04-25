
#ifndef INTARNA_NUSSINOV_HANDLER
#define INTARNA_NUSSINOV_HANDLER

#include "general.h"
#include "RnaSequence.h"

#include <boost/numeric/ublas/triangular.hpp>

namespace IntaRNA {

class NussinovHandler {
public:

  typedef double P_type;  // Probability type
  typedef boost::numeric::ublas::triangular_matrix<P_type, boost::numeric::ublas::upper> P2dMatrix;  // Probability matrix
  typedef boost::numeric::ublas::triangular_matrix<E_type, boost::numeric::ublas::upper> E2dMatrix;  // Energy matrix

  /***
   * Get the partition function Q between the indices (from, to)
   * @param from The start index of the region
   * @param to The end index of the region
   * @param seq The RNA Sequence
   * @param basePairWeight The Boltzmann weight for energy distribution
   * @param minLoopLength The minimum length of loops
   * @param Q Lookup table for Q
   * @param Qb Lookup table for Qb
   * @returns the partition function value
   */
  static E_type getQ(const size_t from, const size_t to, const RnaSequence &seq,
      const E_type basePairWeight, const size_t minLoopLength,
      E2dMatrix &Q, E2dMatrix &Qb);

  /***
   * Get the base partition function Qb between the indices (from, to)
   * @param from The start index of the region
   * @param to The end index of the region
   * @param seq The RNA Sequence
   * @param basePairWeight The Boltzmann weight for energy distribution
   * @param minLoopLength The minimum length of loops
   * @param Q Lookup table for Q
   * @param Qb Lookup table for Qb
   * @returns the partition function value
   */
  static E_type getQb(const size_t from, const size_t to, const RnaSequence &seq,
      const E_type basePairWeight, const size_t minLoopLength,
      E2dMatrix &Q, E2dMatrix &Qb);

  /***
   * Get the base-pair probbility between the indices (from, to)
   * @param from The start index of the region
   * @param to The end index of the region
   * @param seq The RNA Sequence
   * @param basePairWeight The Boltzmann weight for energy distribution
   * @param minLoopLength The minimum length of loops
   * @param Q Lookup table for Q
   * @param Qb Lookup table for Qb
   * @param Pbp Lookup table for base-pair probabilities
   * @returns the probability value
   */
  static P_type getPbp(const size_t from, const size_t to, const RnaSequence &seq,
      const E_type basePairWeight, const size_t minLoopLength,
      E2dMatrix &Q, E2dMatrix &Qb, P2dMatrix &Ppb);

  /***
   * Get the unpaired probbility between the indices (from, to)
   * @param from The start index of the region
   * @param to The end index of the region
   * @param seq The RNA Sequence
   * @param basePairWeight The Boltzmann weight for energy distribution
   * @param minLoopLength The minimum length of loops
   * @param Q Lookup table for Q
   * @param Qb Lookup table for Qb
   * @param Pbp Lookup table for base-pair probabilities
   * @param Pu Lookup table for unpaired probabilities
   * @returns the probability value
   */
  static P_type getPu(const size_t from, const size_t to, const RnaSequence &seq,
      const E_type basePairWeight, const size_t minLoopLength,
      E2dMatrix &Q, E2dMatrix &Qb, P2dMatrix &Ppb, P2dMatrix &Pu);

};

}  // namespace IntaRNA
#endif /* INTARNA_NUSSINOV_HANDLER */

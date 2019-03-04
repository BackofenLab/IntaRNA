
#ifndef INTARNA_NUSSINOV_HANDLER
#define INTARNA_NUSSINOV_HANDLER

#include "IntaRNA/general.h"
#include "IntaRNA/RnaSequence.h"
#include "IntaRNA/Interaction.h"
#include <vector>
#include <utility>

#include <boost/numeric/ublas/triangular.hpp>

namespace IntaRNA {

class NussinovHandler {
public:

  //! Probability triangular matrix
  typedef boost::numeric::ublas::triangular_matrix<Z_type, boost::numeric::ublas::upper> Z2dMatrix;

  //! Energy triangular matrix
  typedef boost::numeric::ublas::triangular_matrix<E_type, boost::numeric::ublas::upper> E2dMatrix;

  //! Index triangular matrix
  typedef boost::numeric::ublas::triangular_matrix<size_t, boost::numeric::ublas::upper> IdxMatrix;

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
  static Z_type getQ(const size_t from, const size_t to, const RnaSequence &seq,
      const Z_type basePairWeight, const size_t minLoopLength,
      Z2dMatrix &Q, Z2dMatrix &Qb);

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
  static Z_type getQb(const size_t from, const size_t to, const RnaSequence &seq,
      const Z_type basePairWeight, const size_t minLoopLength,
      Z2dMatrix &Q, Z2dMatrix &Qb);

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
  static Z_type getPbp(const size_t from, const size_t to, const RnaSequence &seq,
      const Z_type basePairWeight, const size_t minLoopLength,
      Z2dMatrix &Q, Z2dMatrix &Qb, Z2dMatrix &Ppb);

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
  static Z_type getPu(const size_t from, const size_t to, const RnaSequence &seq,
      const Z_type basePairWeight, const size_t minLoopLength,
      Z2dMatrix &Q, Z2dMatrix &Qb, Z2dMatrix &Ppb, Z2dMatrix &Pu);


  /***
   * Get the dotBracket corresponding to the nussinov of the subsequence (from, to)
   * @param from The start of the subsequence
   * @param to The end of the subsequence
   * @param seq The given RNA sequence
   * @param minLoopLength The minimum length of loops
   */
  static std::string dotBracket(const size_t from, const size_t to,
      const RnaSequence &seq, const size_t minLoopLength, const E_type basePairEnergy = Ekcal_2_E(-1.0));

  /***
   * Store all the basepairs in pairs, given a traceback.
   * @param from The start of the pairs
   * @param to The end of the pairs
   * @param traceback The matrix of the traceback indices
   * @param pairs The resulting base-pairs
   */
  static void getBasePairs(const size_t from, const size_t to,
      const IdxMatrix &traceback, Interaction::PairingVec &pairs);

  /***
   * Prints the given matrix to stream. Mainly for debug.
   * @param out the stream to write to
   * @param M the matrix to print
   */
  static void printMatrix( std::ostream & out, const NussinovHandler::E2dMatrix &M);
};

}  // namespace IntaRNA
#endif /* INTARNA_NUSSINOV_HANDLER */

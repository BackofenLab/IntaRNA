#include "AccessibilityBasePair.h"

#include <stdexcept>


namespace IntaRNA {

/////////////////////////////////////////////////////////////////////////////

AccessibilityBasePair::AccessibilityBasePair(const RnaSequence& seq,
    const size_t maxLength, const AccessibilityConstraint * const accConstr_,
    const E_type bpEnergy, const E_type _RT, const size_t minLoopLen) :
      Accessibility(seq, maxLength, accConstr_),
      logPu(seq.size(), seq.size()),
      basePairEnergy(bpEnergy),
      RT(_RT),
      basePairWeight(std::exp(-bpEnergy / _RT)),
      minLoopLength(minLoopLen)
{
  const size_t N = seq.size();
  // create temporary matrices for ED computation
  NussinovHandler::E2dMatrix Q(N, N);
  NussinovHandler::E2dMatrix Qb(N, N);
  NussinovHandler::P2dMatrix Ppb(N, N);
  NussinovHandler::P2dMatrix Pu(N, N);

  logPu.resize(N, N);

  // init temporary matrices
  for (size_t i = 0u; i < N; ++i) {
    for (size_t j = i; j < N; ++j) {
      Q(i, j) = -1.0;
      Qb(i, j) = -1.0;
      Ppb(i, j) = -1.0;
      Pu(i, j) = -1.0;
    }
  }
  // compute ED values
  for (size_t i = 0u; i < N; ++i) {
    for (size_t j = i; j < N; ++j) {
      logPu(i, j) = -RT * std::log(NussinovHandler::getPu(i, j, seq, basePairWeight, minLoopLength, Q, Qb, Ppb, Pu));
    }
  }
}

/////////////////////////////////////////////////////////////////////////////


AccessibilityBasePair::~AccessibilityBasePair()
{}

/////////////////////////////////////////////////////////////////////////////


E_type
AccessibilityBasePair::getED( const size_t from, const size_t to ) const
{
  if (from > to || to < 0 || from >= seq.size()) {
    throw std::runtime_error( "Arguments must satisfy 0 <= from <= to < seq.length" );
  }
  return logPu(from, to);
};

}  // namespace IntaRNA

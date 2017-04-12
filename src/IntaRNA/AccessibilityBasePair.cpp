#include "AccessibilityBasePair.h"

#include <boost/numeric/ublas/triangular.hpp>
#include <stdexcept>


namespace IntaRNA {

/////////////////////////////////////////////////////////////////////////////

AccessibilityBasePair::AccessibilityBasePair(const RnaSequence& seq,
    const size_t maxLength, const AccessibilityConstraint * const accConstr_,
    const E_type bpEnergy, const E_type _RT) :
      Accessibility(seq, maxLength, accConstr_),
      logPu(seq.size(), seq.size()),
      basePairEnergy(bpEnergy),
      RT(_RT)
{
  const size_t N = seq.size();
  E2dMatrix Q(N, N);
  E2dMatrix Qb(N, N);
  P2dMatrix Ppb(N, N);
  P2dMatrix Pu(N, N);

  logPu.resize(N, N);

  for (size_t i = 0u; i < N; ++i) {
    for (size_t j = i; j < N; ++j) {
      Q(i, j) = -1.0;
      Qb(i, j) = -1.0;
      Ppb(i, j) = -1.0;
      Pu(i, j) = -1.0;
    }
  }
  for (size_t i = 0u; i < N; ++i) {
    for (size_t j = i; j < N; ++j) {
      logPu(i, j) = -RT * std::log(getPu(i, j, Q, Qb, Ppb, Pu));
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


/////////////////////////////////////////////////////////////////////////////


E_type
AccessibilityBasePair::getQ(const size_t i, const size_t j,
    AccessibilityBasePair::E2dMatrix &Q, AccessibilityBasePair::E2dMatrix &Qb) {
  if (i > j || i < 0 || j >= seq.size() || i + minLoopLength >= j) {
    return 1.0;
  }
  E_type &ret = Q(i, j);
  // If value is already computed, return it
  if (ret > -0.5) {
    return ret;
  }
  // Else compute Q
  ret = getQ(i, j - 1, Q, Qb);
  for (size_t k = i; k + minLoopLength < j; ++k) {
      ret += getQ(i, k - 1, Q, Qb) * getQb(k, j, Q, Qb);
  }
  return ret;
}

/////////////////////////////////////////////////////////////////////////////


E_type
AccessibilityBasePair::getQb(const size_t i, const size_t j,
    AccessibilityBasePair::E2dMatrix &Q, AccessibilityBasePair::E2dMatrix &Qb) {
  if (i > j || i < 0 || j >= seq.size() || i + minLoopLength >= j) {
    return 0.0;
  }
  E_type &ret = Qb(i, j);
  // If value is already computed, return it
  if (ret > -0.5) {
    return ret;
  }
  // Else compute Qb
  if (RnaSequence::areComplementary(seq, seq, i, j)) {
    ret = getQ(i + 1, j - 1, Q, Qb) * basePairWeight;
  } else {
    ret = 0;
  }
  return ret;
}

/////////////////////////////////////////////////////////////////////////////


AccessibilityBasePair::P_type
AccessibilityBasePair::getPbp(const size_t i, const size_t j,
    AccessibilityBasePair::E2dMatrix &Q, AccessibilityBasePair::E2dMatrix &Qb,
    AccessibilityBasePair::P2dMatrix &Ppb) {
  if (i > j || i < 0 || j >= seq.size() || i + minLoopLength >= j) {
    return 0.0;
  }
  P_type &ret = Ppb(i, j);
  // If value is already computed, return it
  if (ret > -0.5) {
    return ret;
  }
  // Else compute Pbp
  ret = (getQb(i, j, Q, Qb) * getQ(0, i - 1, Q, Qb) *
         getQ(j + 1, seq.size() - 1, Q, Qb) / getQ(0, seq.size() - 1, Q, Qb));
  for (size_t p = 0; p < i; ++p) {
    for (size_t q = j + 1; q < seq.size(); ++q) {
      if (RnaSequence::areComplementary(seq, seq, p, q)) {
        ret += (basePairWeight * getPbp(p, q, Q, Qb, Ppb) *
                getQ(p + 1, i - 1, Q, Qb) * getQb(i, j, Q, Qb) *
                getQ(j + 1, q - 1, Q, Qb) / getQb(p, q, Q, Qb));
      }
    }
  }
  return ret;
}

/////////////////////////////////////////////////////////////////////////////


AccessibilityBasePair::P_type
AccessibilityBasePair::getPu(const size_t i, const size_t j,
    AccessibilityBasePair::E2dMatrix &Q, AccessibilityBasePair::E2dMatrix &Qb,
    AccessibilityBasePair::P2dMatrix &Pbp, AccessibilityBasePair::P2dMatrix &Pu) {
  if (i > j || i < 0 || j >= seq.size() || i + minLoopLength >= j) {
    return 0.0;
  }
  P_type &ret = Pu(i, j);
  // If value is already computed, return it
  if (ret > -0.5) {
    return ret;
  }
  // Else compute Pu
  ret = getQ(0, i - 1, Q, Qb) * getQ(j + 1, seq.size() - 1, Q, Qb) / getQ(0, seq.size() - 1, Q, Qb);
  for (size_t p = 0; p < i; ++p) {
    for (size_t q = j + 1; q < seq.size(); ++q) {
      if (RnaSequence::areComplementary(seq, seq, p, q)) {
        ret += (basePairWeight * getPbp(p, q, Q, Qb, Pbp) *
                getQ(p + 1, i - 1, Q, Qb) * getQ(j + 1, q - 1, Q, Qb) /
                getQb(p, q, Q, Qb));
      }
    }
  }
  return ret;
}

}  // namespace IntaRNA

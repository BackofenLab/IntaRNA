
#include "IntaRNA/InteractionEnergyBasePair.h"

namespace IntaRNA {

void InteractionEnergyBasePair::computeES(const Accessibility &accS,
    InteractionEnergyBasePair::E2dMatrix &logQ) {
  const size_t N = accS.getSequence().size();

  E2dMatrix Q(N, N);
  E2dMatrix Qb(N, N);

  logQ.resize(N, N);
  for (size_t i = 0u; i < N; ++i) {
    for (size_t j = 0u; j < N; ++j) {
      Q(i, j) = -1.0;
      Qb(i, j) = -1.0;
    }
  }

  for (size_t i = 0u; i < N; ++i) {
    for (size_t j = i; j < N; ++j) {
      logQ(i, j) = -RT * std::log(getQ(i, j, accS.getSequence(), Q, Qb));
    }
  }
}


E_type
InteractionEnergyBasePair::getQ(const size_t i, const size_t j, const RnaSequence &seq,
    InteractionEnergyBasePair::E2dMatrix &Q, InteractionEnergyBasePair::E2dMatrix &Qb) {
  if (i > j || i < 0 || j >= seq.size() || i + minLoopLength >= j) {
    return 1.0;
  }
  E_type &ret = Q(i, j);
  if (ret > -0.5) {
    return ret;
  }
  ret = getQ(i, j - 1, seq, Q, Qb);
  for (size_t k = i; k + minLoopLength < j; ++k) {
      ret += getQ(i, k - 1, seq, Q, Qb) * getQb(k, j, seq, Q, Qb);
  }
  return ret;
}


E_type
InteractionEnergyBasePair::getQb(const size_t i, const size_t j, const RnaSequence &seq,
    InteractionEnergyBasePair::E2dMatrix &Q, InteractionEnergyBasePair::E2dMatrix &Qb) {
  if (i > j || i < 0 || j >= seq.size() || i + minLoopLength >= j) {
    return 0.0;
  }
  E_type &ret = Qb(i, j);
  if (ret > -0.5) {
    return ret;
  }
  if (RnaSequence::areComplementary(seq, seq, i, j)) {
    ret = getQ(i + 1, j - 1, seq, Q, Qb) * boltzmann;
  } else {
    ret = 0;
  }
  return ret;
}

} // namespace IntaRNA

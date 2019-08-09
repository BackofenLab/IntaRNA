
#include "IntaRNA/InteractionEnergyBasePair.h"

namespace IntaRNA {

void InteractionEnergyBasePair::computeES(const RnaSequence &seq,
    NussinovHandler::E2dMatrix &logQ) {
  const size_t N = seq.size();

  NussinovHandler::Z2dMatrix Q(N, N);
  NussinovHandler::Z2dMatrix Qb(N, N);

  logQ.resize(N, N);
  for (size_t i = 0u; i < N; ++i) {
    for (size_t j = i; j < N; ++j) {
      Q(i, j) = -1.0;
      Qb(i, j) = -1.0;
    }
  }

  for (size_t i = 0u; i < N; ++i) {
    for (size_t j = i; j < N; ++j) {
      const Z_type q_val = NussinovHandler::getQ(i, j, seq, basePairWeight, minLoopLength, Q, Qb);
      if (Z_equal(q_val, 1.0)) {
        logQ(i, j) = E_INF;
      } else {
        logQ(i, j) = getE(q_val);
      }
    }
  }
}

}  // namespace IntaRNA


#include "NussinovHandler.h"

namespace IntaRNA {

E_type
NussinovHandler::getQ(const size_t i, const size_t j, const RnaSequence &seq,
    const E_type bpWeight, const size_t minLoopLength,
    NussinovHandler::E2dMatrix &Q, NussinovHandler::E2dMatrix &Qb) {
  if (i > j || i < 0 || j >= seq.size() || i + minLoopLength >= j) {
    return 1.0;
  }
  E_type &ret = Q(i, j);
  // If value is already computed, return it
  if (ret > -0.5) {
    return ret;
  }
  // Else compute Q
  ret = getQ(i, j - 1, seq, bpWeight, minLoopLength, Q, Qb);
  for (size_t k = i; k + minLoopLength < j; ++k) {
      ret += (getQ(i, k - 1, seq, bpWeight, minLoopLength, Q, Qb) *
              getQb(k, j, seq, bpWeight, minLoopLength, Q, Qb));
  }
  return ret;
}


E_type
NussinovHandler::getQb(const size_t i, const size_t j, const RnaSequence &seq,
    const E_type bpWeight, const size_t minLoopLength,
    NussinovHandler::E2dMatrix &Q, NussinovHandler::E2dMatrix &Qb) {
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
    ret = getQ(i + 1, j - 1, seq, bpWeight, minLoopLength, Q, Qb) * bpWeight;
  } else {
    ret = 0;
  }
  return ret;
}


NussinovHandler::P_type
NussinovHandler::getPbp(const size_t i, const size_t j, const RnaSequence &seq,
    const E_type bpWeight, const size_t minLoopLength,
    NussinovHandler::E2dMatrix &Q, NussinovHandler::E2dMatrix &Qb,
    NussinovHandler::P2dMatrix &Ppb) {
  if (i > j || i < 0 || j >= seq.size() || i + minLoopLength >= j) {
    return 0.0;
  }
  P_type &ret = Ppb(i, j);
  // If value is already computed, return it
  if (ret > -0.5) {
    return ret;
  }
  // Else compute Pbp
  ret = (getQb(i, j, seq, bpWeight, minLoopLength, Q, Qb) *
         getQ(0, i - 1, seq, bpWeight, minLoopLength, Q, Qb) *
         getQ(j + 1, seq.size() - 1, seq, bpWeight, minLoopLength, Q, Qb) /
         getQ(0, seq.size() - 1, seq, bpWeight, minLoopLength, Q, Qb));
  for (size_t p = 0; p < i; ++p) {
    for (size_t q = j + 1; q < seq.size(); ++q) {
      if (RnaSequence::areComplementary(seq, seq, p, q)) {
        ret += (bpWeight *
                getPbp(p, q, seq, bpWeight, minLoopLength, Q, Qb, Ppb) *
                getQ(p + 1, i - 1, seq, bpWeight, minLoopLength, Q, Qb) *
                getQb(i, j, seq, bpWeight, minLoopLength, Q, Qb) *
                getQ(j + 1, q - 1, seq, bpWeight, minLoopLength, Q, Qb) /
                getQb(p, q, seq, bpWeight, minLoopLength, Q, Qb));
      }
    }
  }
  return ret;
}


NussinovHandler::P_type
NussinovHandler::getPu(const size_t i, const size_t j, const RnaSequence &seq,
    const E_type bpWeight, const size_t minLoopLength,
    NussinovHandler::E2dMatrix &Q, NussinovHandler::E2dMatrix &Qb,
    NussinovHandler::P2dMatrix &Pbp, NussinovHandler::P2dMatrix &Pu) {
  if (i > j || i < 0 || j >= seq.size() || i + minLoopLength >= j) {
    return 0.0;
  }
  P_type &ret = Pu(i, j);
  // If value is already computed, return it
  if (ret > -0.5) {
    return ret;
  }
  // Else compute Pu
  ret = (getQ(0, i - 1, seq, bpWeight, minLoopLength, Q, Qb) *
         getQ(j + 1, seq.size() - 1, seq, bpWeight, minLoopLength, Q, Qb) /
         getQ(0, seq.size() - 1, seq, bpWeight, minLoopLength, Q, Qb));
  for (size_t p = 0; p < i; ++p) {
    for (size_t q = j + 1; q < seq.size(); ++q) {
      if (RnaSequence::areComplementary(seq, seq, p, q)) {
        ret += (bpWeight *
                getPbp(p, q, seq, bpWeight, minLoopLength, Q, Qb, Pbp) *
                getQ(p + 1, i - 1, seq, bpWeight, minLoopLength, Q, Qb) *
                getQ(j + 1, q - 1, seq, bpWeight, minLoopLength, Q, Qb) /
                getQb(p, q, seq, bpWeight, minLoopLength, Q, Qb));
      }
    }
  }
  return ret;
}


}  // namespace IntaRNA

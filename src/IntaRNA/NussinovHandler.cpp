
#include "IntaRNA/NussinovHandler.h"

namespace IntaRNA {

void
NussinovHandler::getBasePairs(
    const size_t from, 
    const size_t to,
    const IdxMatrix &traceback, 
    std::vector<std::pair<size_t, size_t>> &pairs) 
{
  if (from >= to) {
    return ;
  }
  if (traceback(from, to) == to) {
    getBasePairs(from, to - 1, traceback, pairs);
  } else {
    size_t k = traceback(from, to);
    if (from + 1 <= k)  {
      getBasePairs(from, k - 1, traceback, pairs);
    }
    if (k + 2 <= to) {
      getBasePairs(k + 1, to - 1, traceback, pairs);
    }

    pairs.push_back(std::pair<size_t, size_t>(k, to));
  }
}

std::string
NussinovHandler::dotBracket(const size_t from, const size_t to,
    const RnaSequence &seq, const size_t minLoopLength) 
{
  const size_t len = to - from + 1, offset = from;
  NussinovHandler::E2dMatrix nussinov(len + 1, len + 1);
  NussinovHandler::IdxMatrix traceback(len + 1, len + 1);
  for (size_t i = 0; i < len; ++i)
    for (size_t j = i; j < len; ++j)
      nussinov(i, j) = 0, traceback(i, j) = j;

  for (size_t dis = 1; dis < len; ++dis) {
    for (size_t beg = 0; beg + dis < len; ++beg) {
      const size_t i = beg, j = beg + dis;
      nussinov(i, j) = nussinov(i, j - 1);
      traceback(i, j) = j;
      for (size_t k = i; k + minLoopLength < j; ++k) {
        if (RnaSequence::areComplementary(seq, seq, offset + k, offset + j)) {
          const E_type val = ((i + 1 <= k)  ? nussinov(i , k - 1) : 0) + ((k + 2 <= j) ? nussinov(k + 1, j - 1) : 0) + 1;
          if (val > nussinov(i, j)) {
            nussinov(i, j) = val;
            traceback(i, j) = k;
          }
        }
      }
    }
  }
  std::string result(len, '.');
  std::vector<std::pair<size_t, size_t>> basepairs;
  getBasePairs(0, len - 1, traceback, basepairs);
  for (size_t k = 0; k < basepairs.size(); ++k) {
    size_t i = basepairs[k].first, j = basepairs[k].second;
    // assert(result[i] == '.' && result[j] == '.');
    result[i] = '(';
    result[j] = ')';
  }
  return result;
}

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

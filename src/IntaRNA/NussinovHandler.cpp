
#include "IntaRNA/NussinovHandler.h"

namespace IntaRNA {

void
NussinovHandler::getBasePairs(
    const size_t from,
    const size_t to,
    const IdxMatrix &traceback,
    Interaction::PairingVec &pairs)
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
    const RnaSequence &seq, const size_t minLoopLength, const E_type basePairEnergy)
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
          const E_type val = ((i + 1 <= k)  ? nussinov(i , k - 1) : 0) + ((k + 2 <= j) ? nussinov(k + 1, j - 1) : 0) + basePairEnergy;
          if (val < nussinov(i, j)) {
            nussinov(i, j) = val;
            traceback(i, j) = k;
          }
        }
      }
    }
  }
  std::string result(len, '.');
  Interaction::PairingVec basepairs;
  getBasePairs(0, len - 1, traceback, basepairs);
  for (size_t k = 0; k < basepairs.size(); ++k) {
    size_t i = basepairs[k].first, j = basepairs[k].second;
    // assert(result[i] == '.' && result[j] == '.');
    result[i] = '(';
    result[j] = ')';
  }
  return result;
}

Z_type
NussinovHandler::getQ(const size_t i, const size_t j, const RnaSequence &seq,
    const Z_type bpWeight, const size_t minLoopLength,
    NussinovHandler::Z2dMatrix &Q, NussinovHandler::Z2dMatrix &Qb) {
  if (i > j || j >= seq.size()) {
    return 1.0;
  }
  Z_type &ret = Q(i, j);
  // If value is already computed, return it
  if (ret > -0.5) {
    return ret;
  }
  // Else compute Q
  ret = (j>0?getQ(i, j - 1, seq, bpWeight, minLoopLength, Q, Qb):1.0);
  for (size_t k = i; k + minLoopLength < j; ++k) {
      ret += ( (k>0?getQ(i, k - 1, seq, bpWeight, minLoopLength, Q, Qb):1.0) *
              getQb(k, j, seq, bpWeight, minLoopLength, Q, Qb));
  }
  return ret;
}


Z_type
NussinovHandler::getQb(const size_t i, const size_t j, const RnaSequence &seq,
    const Z_type bpWeight, const size_t minLoopLength,
    NussinovHandler::Z2dMatrix &Q, NussinovHandler::Z2dMatrix &Qb) {
  if (j >= seq.size()) {
    return 1.0;
  }
  if (i + minLoopLength >= j) {
    return 0.0;
  }
  Z_type &ret = Qb(i, j);
  // If value is already computed, return it
  if (ret > -0.5) {
    return ret;
  }
  // Else compute Qb
  if (RnaSequence::areComplementary(seq, seq, i, j)) {
    ret = getQ(i + 1, j - 1, seq, bpWeight, minLoopLength, Q, Qb) * bpWeight;
  } else {
    ret = 0.0;
  }
  return ret;
}


Z_type
NussinovHandler::getPbp(const size_t i, const size_t j, const RnaSequence &seq,
    const Z_type bpWeight, const size_t minLoopLength,
    NussinovHandler::Z2dMatrix &Q, NussinovHandler::Z2dMatrix &Qb,
    NussinovHandler::Z2dMatrix &Ppb) {
  if (j >= seq.size() || i + minLoopLength >= j) {
    return 0.0;
  }
  Z_type &ret = Ppb(i, j);
  // If value is already computed, return it
  if (ret > -0.5) {
    return ret;
  }
  // Else compute Pbp
  ret = (getQb(i, j, seq, bpWeight, minLoopLength, Q, Qb) *
         (i>0?getQ(0, i - 1, seq, bpWeight, minLoopLength, Q, Qb):1.0) *
         getQ(j + 1, seq.size() - 1, seq, bpWeight, minLoopLength, Q, Qb) /
         getQ(0, seq.size() - 1, seq, bpWeight, minLoopLength, Q, Qb));
  for (size_t p = 0; p < i; ++p) {
    for (size_t q = j + 1; q < seq.size(); ++q) {
      if (p+minLoopLength < q && RnaSequence::areComplementary(seq, seq, p, q)) {
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


Z_type
NussinovHandler::getPu(const size_t i, const size_t j, const RnaSequence &seq,
    const Z_type bpWeight, const size_t minLoopLength,
    NussinovHandler::Z2dMatrix &Q, NussinovHandler::Z2dMatrix &Qb,
    NussinovHandler::Z2dMatrix &Pbp, NussinovHandler::Z2dMatrix &Pu) {
  if (i > j || j >= seq.size()) {
    return 0.0;
  }
  Z_type &ret = Pu(i, j);
  // If value is already computed, return it
  if (ret > -0.5) {
    return ret;
  }
  // Else compute Pu
  ret = ((i>0?getQ(0, i - 1, seq, bpWeight, minLoopLength, Q, Qb):1.0) *
         getQ(j + 1, seq.size() - 1, seq, bpWeight, minLoopLength, Q, Qb) /
         getQ(0, seq.size() - 1, seq, bpWeight, minLoopLength, Q, Qb));
  for (size_t p = 0; p < i; ++p) {
    for (size_t q = j + 1; q < seq.size(); ++q) {
      if (p+minLoopLength<q && RnaSequence::areComplementary(seq, seq, p, q)) {
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

void
NussinovHandler::
printMatrix( std::ostream & out, const NussinovHandler::E2dMatrix &M)
{
    std::cout <<"\n\n########\n";
    for (int i =0; i<M.size1(); i++) {
    	std::cout <<(i);
    	for (int j=0; j <M.size2(); j++) {
    		if (j<i) {
    			std::cout <<" ---";
    		} else {
    			std::cout <<" "<<M(i,j);
    		}
    	}
    	std::cout <<std::endl;
    }
    std::cout <<"########\n";
}


}  // namespace IntaRNA

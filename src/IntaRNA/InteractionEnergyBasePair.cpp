
#include "IntaRNA/InteractionEnergyBasePair.h"

namespace IntaRNA {

////////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyBasePair::
getEall1() const
{
	// compute Z if needed
	if (E_isINF(Eall1)) {
		Eall1 = computeIntraEall( accS1 );
	}
	return Eall1;
}

////////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyBasePair::
getEall2() const
{
	// compute Z if needed
	if (E_isINF(Eall2)) {
		Eall2 = computeIntraEall( accS2.getAccessibilityOrigin() );
	}
	return Eall2;
}

////////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyBasePair::
computeIntraEall( const Accessibility & acc ) const
{
	if ( !acc.getAccConstraint().isEmpty() ) {
		INTARNA_NOT_IMPLEMENTED("InteractionEnergyBasePair: accessibility constraints for ensemble energy computation not supported");
	}

	const size_t N = acc.getSequence().size();

	// check if any base pair can be formed
	if (N < minLoopLength) {
		return 0;
	}

	// create temporary matrices for ED computation
	NussinovHandler::Z2dMatrix Q(N, N);
	NussinovHandler::Z2dMatrix Qb(N, N);

	// init temporary matrices
	for (size_t i = 0u; i < N; ++i) {
		for (size_t j = i; j < N; ++j) {
			Q(i, j) = -1.0;
			Qb(i, j) = -1.0;
		}
	}

	// compute partition function and convert to ensemble energy
	return getE( NussinovHandler::getQ( 0, N-1, acc.getSequence(), basePairWeight, minLoopLength, Q, Qb) );

}

////////////////////////////////////////////////////////////////////////////

}  // namespace IntaRNA

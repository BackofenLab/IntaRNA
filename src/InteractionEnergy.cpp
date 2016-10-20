
#include "InteractionEnergy.h"

#include <cmath>

#include <iostream>

////////////////////////////////////////////////////////////////////////////

InteractionEnergy::InteractionEnergy( const Accessibility & accS1
				, const ReverseAccessibility & accS2
				, const size_t maxInternalLoopSize1
				, const size_t maxInternalLoopSize2
		)
  :
	accS1(accS1)
	, accS2(accS2)
	, maxInternalLoopSize1(maxInternalLoopSize1)
	, maxInternalLoopSize2(maxInternalLoopSize2)

{
}

////////////////////////////////////////////////////////////////////////////

InteractionEnergy::~InteractionEnergy()
{
}

////////////////////////////////////////////////////////////////////////////

bool
InteractionEnergy::
isAllowedLoopRegion( const RnaSequence& seq, const size_t i, const size_t j, const size_t maxInternalLoopSize )
{
	// ensure index and loop size validity
	return	   i < seq.size()
			&& j < seq.size()
			&& seq.asString().at(i) != 'N'
			&& seq.asString().at(j) != 'N'
			&& i <= j
			&& (j-i) <= (1+maxInternalLoopSize);

}

////////////////////////////////////////////////////////////////////////////

bool
InteractionEnergy::
areComplementary( const size_t i1, const size_t i2 ) const
{
	return RnaSequence::areComplementary( accS1.getSequence(), accS2.getSequence(), i1, i2);
}

////////////////////////////////////////////////////////////////////////////

bool
InteractionEnergy::
isValidInternalLoop( const size_t i1, const size_t j1, const size_t i2, const size_t j2 ) const
{
	return
		   (j1-i1>0 && j2-i2>0)
		&& areComplementary( i1, i2)
		&& areComplementary( j1, j2)
		&& InteractionEnergy::isAllowedLoopRegion(accS1.getSequence(), i1, j1, maxInternalLoopSize1)
		&& InteractionEnergy::isAllowedLoopRegion(accS2.getSequence(), i2, j2, maxInternalLoopSize2)
		;
}

////////////////////////////////////////////////////////////////////////////

const Accessibility &
InteractionEnergy::
getAccessibility1() const
{
	return accS1;
}

////////////////////////////////////////////////////////////////////////////

const ReverseAccessibility &
InteractionEnergy::
getAccessibility2() const
{
	return accS2;
}

////////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergy::
getPr_danglingLeft( const size_t i1, const size_t j1, const size_t i2, const size_t j2 ) const
{
	// initial probabilities
	E_type probDangle1 = 1.0, probDangle2 = 1.0;

	// if dangle1 possible
	if (i1>0)  {
		// Pr( i1-1 is unpaired | i1..j1 unpaired )
		probDangle1 =
			std::max( (E_type)0.0
					, std::min( (E_type)1.0
							, getBoltzmannWeight(
									getAccessibility1().getED( i1-1, j1 )
									- getAccessibility1().getED( i1, j1 ) )
							)
					)
			;
	}
	// if dangle2 possible
	if (i2>0)  {
		// Pr( i2-1 is unpaired | i2..j2 unpaired )
		probDangle2 =
			std::max( (E_type)0.0
					, std::min( (E_type)1.0
							, getBoltzmannWeight(
									getAccessibility2().getED( i2-1, j2 )
									- getAccessibility2().getED( i2, j2 ) )
							)
					)
			;
	}

	// get overall probability
	return probDangle1 * probDangle2;
}

////////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergy::
getPr_danglingRight( const size_t i1, const size_t j1, const size_t i2, const size_t j2 ) const
{
	// initial probabilities
	E_type probDangle1 = 1.0, probDangle2 = 1.0;

	// if dangle1 possible
	if (j1+1<getAccessibility1().getSequence().size())  {
		// Pr( j1+1 is unpaired | i1..j1 unpaired )
		probDangle1 =
			std::max( (E_type)0.0
					, std::min( (E_type)1.0
							, getBoltzmannWeight(
									getAccessibility1().getED( i1, j1+1 )
									- getAccessibility1().getED( i1, j1 ) )
							)
					)
			;
	}
	// if dangle2 possible
	if (j2+1<getAccessibility2().getSequence().size())  {
		// Pr( j2+1 is unpaired | i2..j2 unpaired )
		probDangle2 =
			std::max( (E_type)0.0
					, std::min( (E_type)1.0
							, getBoltzmannWeight(
									getAccessibility2().getED( i2, j2+1 )
									- getAccessibility2().getED( i2, j2 ) )
							)
					)
			;
	}

	// get overall probability
	return probDangle1 * probDangle2;
}

////////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergy::
getE( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2
		, const E_type hybridE ) const
{
	// check if hybridization energy is not infinite
	if ( E_isNotINF(hybridE) ) {
		// compute overall interaction energy
		return hybridE
				// accessibility penalty
				+ getAccessibility1().getED( i1, j1 )
				+ getAccessibility2().getED( i2, j2 )
				// dangling end penalty
				// weighted by the probability that ends are unpaired
				+ (getE_danglingLeft( i1, i2 )*getPr_danglingLeft(i1,j1,i2,j2))
				+ (getE_danglingRight( j1, j2 )*getPr_danglingRight(i1,j1,i2,j2))
				// helix closure penalty
				+ getE_endLeft( i1, i2 )
				+ getE_endRight( j1, j2 )
				;
	} else {
		// hybridE is infinite, thus overall energy is infinity as well
		return E_INF;
	}
}

////////////////////////////////////////////////////////////////////////////

InteractionEnergy::
EnergyContributions
InteractionEnergy::
getE_contributions( const Interaction & interaction ) const
{

	// temporary access to range indices
	const size_t i1 = interaction.basePairs.begin()->first;
	const size_t i2 = getAccessibility2().getReversedIndex(interaction.basePairs.begin()->second);
	const size_t j1 = interaction.basePairs.rbegin()->first;
	const size_t j2 = getAccessibility2().getReversedIndex(interaction.basePairs.rbegin()->second);

	// fill contribution data structure
	EnergyContributions contr;
	contr.init = getE_init();
	contr.ED1 = getAccessibility1().getED( i1, j1 );
	contr.ED2 = getAccessibility2().getED( i2, j2 );
	contr.dangleLeft = (getE_danglingLeft( i1, i2 )*getPr_danglingLeft(i1,j1,i2,j2));
	contr.dangleRight = (getE_danglingRight( j1, j2 )*getPr_danglingRight(i1,j1,i2,j2));
	contr.endLeft = getE_endLeft( i1, i2 );
	contr.endRight = getE_endRight( j1, j2 );
	// compute loop energy
	contr.loops = interaction.energy
					- contr.init
					- contr.ED1
					- contr.ED2
					- contr.dangleLeft
					- contr.dangleRight
					- contr.endLeft
					- contr.endRight
					;

	// final contribution distribution
	return contr;
}

////////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergy::
getBoltzmannWeight( const E_type e ) const
{
	// TODO can be optimized when using exp-energies from VRNA
	return std::exp( - e / getRT() );
}

////////////////////////////////////////////////////////////////////////////



#include "InteractionEnergy.h"

#include <cmath>

#include <iostream>

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
							, getBoltzmannWeight( getED1(i1-1,j1)-getED1(i1,j1) )
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
							, getBoltzmannWeight( getED2(i2-1,j2)-getED2(i2,j2) )
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
	if (j1+1<size1())  {
		// Pr( j1+1 is unpaired | i1..j1 unpaired )
		probDangle1 =
			std::max( (E_type)0.0
					, std::min( (E_type)1.0
							, getBoltzmannWeight( getED1(i1,j1+1)-getED1(i1,j1) )
							)
					)
			;
	}
	// if dangle2 possible
	if (j2+1<size2())  {
		// Pr( j2+1 is unpaired | i2..j2 unpaired )
		probDangle2 =
			std::max( (E_type)0.0
					, std::min( (E_type)1.0
							, getBoltzmannWeight( getED2(i2,j2+1)-getED2(i2,j2) )
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
				+ getED1( i1, j1 )
				+ getED2( i2, j2 )
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
	contr.ED1 = getED1( i1, j1 );
	contr.ED2 = getED2( i2, j2 );
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


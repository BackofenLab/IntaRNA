
#include "InteractionEnergy.h"

#include <cmath>

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
		   areComplementary( i1, i2)
		&& areComplementary( j1, j2)
		&& InteractionEnergy::isAllowedLoopRegion(accS1.getSequence(), i1, j1, maxInternalLoopSize1)
		&& InteractionEnergy::isAllowedLoopRegion(accS2.getSequence(), i2, j2, maxInternalLoopSize2)
		&& ( (j1-i1==0 && j2-i2==0) || (j1-i1>0 && j2-i2>0) )
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
getE( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2
		, const E_type hybridE ) const
{
	return hybridE
			// accessibility penalty
			+ getAccessibility1().getED( i1, j1 )
			+ getAccessibility2().getED( i2, j2 )
			// dangling end penalty
			// weighted by the probability that ends are unpaired (not included in RNAup)
			+ (getDanglingLeft( i1, i2 )
					// Pr( i1-1 is unpaired | i1..j1 unpaired )
					*getBoltzmannWeight( getAccessibility1().getED( (i1==0?i1:i1-1), j1 ) - getAccessibility1().getED( i1, j1 ) )
					// Pr( i2-1 is unpaired | i2..j2 unpaired )
					*getBoltzmannWeight( getAccessibility2().getED( (i2==0?i2:i2-1), j2 ) - getAccessibility2().getED( i2, j2 ) )
					)
			+ (getDanglingRight( j1, j2 )
					// Pr( j1+1 is unpaired | i1..j1 unpaired )
					*getBoltzmannWeight( getAccessibility1().getED( i1, std::min(getAccessibility1().getSequence().size()-1,j1+1) ) - getAccessibility1().getED( i1, j1 ) )
					// Pr( j2+1 is unpaired | i2..j2 unpaired )
					*getBoltzmannWeight( getAccessibility2().getED( i2, std::min(getAccessibility2().getSequence().size()-1,j2+1) ) - getAccessibility2().getED( i2, j2 ) )
					)
			// helix closure penalty (not included in RNAup)
			+ getEndLeft( i1, i2 )
			+ getEndRight( j1, j2 )
			;
}

////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergy::
getBoltzmannWeight( const E_type e ) const
{
	// TODO can be optimized when using exp-energies from VRNA
	return std::exp( - e / getRT() );
}

////////////////////////////////////////////////////////////////////////////


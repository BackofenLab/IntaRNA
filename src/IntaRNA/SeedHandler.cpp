
#include "IntaRNA/SeedHandler.h"

namespace IntaRNA {


////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

bool
SeedHandler::
isFeasibleSeedBasePair( const size_t i1, const size_t i2, const bool atEndOfSeed ) const
{

	return		i1 < energy.size1() && i2 < energy.size2()
			&&	energy.isAccessible1(i1)
			&&	energy.isAccessible2(i2)
			&&	energy.areComplementary(i1,i2)
			&&	seedConstraint.getMaxED() >= energy.getED1( i1,i1 )
			&&	seedConstraint.getMaxED() >= energy.getED2( i2,i2 )
			&&  (seedConstraint.isGUallowed() || (!energy.isGU( i1, i2 ))) // check for GU bp constraint
			&&  (!atEndOfSeed || (seedConstraint.isGUendAllowed() || (!energy.isGU( i1, i2 )))) // check for GU ends
			&&	(seedConstraint.getRanges1().empty() || seedConstraint.getRanges1().covers(i1))
			&&	(seedConstraint.getRanges2().empty() || seedConstraint.getRanges2().covers(i2))
			;
}

//////////////////////////////////////////////////////////////////////////

bool
SeedHandler::
updateToNextSeed( size_t & i1_out, size_t & i2_out
		, const size_t i1min, const size_t i1max
		, const size_t i2min, const size_t i2max
		) const
{
	size_t i1=i1_out, i2=i2_out;
	// find true max value
	const size_t i1maxVal = std::min(energy.size1()-1,i1max)
				, i2maxVal = std::min(energy.size2()-1,i2max);
	// find first seed if out of bound or no valid seed boundary
	if ( i1 < i1min || i1 > i1maxVal
		|| i2 < i2min || i2 > i2maxVal
		|| !(isSeedBound(i1,i2)) )
	{
		// first potential seed position
		i1 = i1min;
		i2 = i2min;
	} else {
		// update to next potential seed position
		if (++i1 > i1maxVal) {
			i1 = i1min;
			i2++;
		}
	}

	// find next valid seed start within range
	while( i2 <= i2maxVal && !(isSeedBound(i1,i2))) {
		// update seed position within range
		if (++i1 > i1maxVal) {
			i1 = i1min;
			i2++;
		}
	}

	// check if we found a valid seed in the range
	if (i1 <= i1maxVal && i2 <= i2maxVal) {
		i1_out = i1;
		i2_out = i2;
		return true;
	}
	// no valid next seed found
	return false;
}

//////////////////////////////////////////////////////////////////////////

void
SeedHandler::
addSeeds( Interaction & i ) const
{
	// screen all base pairs
	for ( Interaction::PairingVec::const_iterator bp = i.basePairs.begin(); bp != i.basePairs.end(); bp++ ) {
		// get internal indices
		size_t i1 = energy.getIndex1(*bp), i2 = energy.getIndex2(*bp);
		// check if seed boundary
		if (isSeedBound(i1,i2)) {
			// check if seed end is within interaction i
			size_t j1 = i1+getSeedLength1(i1,i2)-1, j2 = i2+getSeedLength2(i1,i2)-1;
			if ( std::find(i.basePairs.begin(), i.basePairs.end(), energy.getBasePair(j1,j2)) == i.basePairs.end() ) {
				// seed end not in i's base pairs -> skip this seed
				continue;
			}
			// seed end is covered -> traceback and check remaining base pairs
			Interaction tmp( *(i.s1), *(i.s2) );
			traceBackSeed( tmp, i1, i2 );
			Interaction::PairingVec::const_iterator bpS(bp);
			for ( auto & sBP : tmp.basePairs ) {
				bpS++;
				// check for seed base pair
				if ( sBP != (*bpS) ) {
					continue;
				}
			}
			// add seed to interaction
			if (i.seed == NULL) { i.seed = new Interaction::SeedSet(); }
			i.seed->insert( Interaction::Seed( *bp, energy.getBasePair(j1,j2)
					// compute overall energy
					, energy.getE(i1,j1,i2,j2,getSeedE(i1,i2))+energy.getE_init()
					) );
		}
	}
}


//////////////////////////////////////////////////////////////////////////

bool
SeedHandler::
areLoopOverlapping( const size_t i1, const size_t i2
				, const size_t k1, const size_t k2 ) const
{
	// check if valid seeds
	if (isSeedBound(i1,i2) && isSeedBound(k1,k2)) {
		// check for identity
		if (i1 == k1 && i2 == k2 ) {
			return true;
		}
		// ensure order x < y
		size_t x1=i1, x2=i2, y1=k1, y2=k2;
		if (i1 > k1) {
			x1=k1; x2=k2; y1=i1; y2=i2;
		}
		// ensure order of seq2 indices; otherwise not overlapping
		if (x2 >= y2) {
			return false;
		}
		// check if ranges are overlapping at all
		if ( ((x1+getSeedLength1(x1,x2)-1) <= y1)
			|| ((x2+getSeedLength2(x1,x2)-1) <= y2) )
		{
			return false;
		}
		// trace seed base pairs of seedX (excluding first)
		Interaction xBP( energy.getAccessibility1().getSequence(), energy.getAccessibility2().getAccessibilityOrigin().getSequence() );
		traceBackSeed( xBP, x1, x2 );
		xBP.basePairs.push_back( energy.getBasePair( x1 + getSeedLength1(x1,x2) - 1, x2 + getSeedLength2(x1,x2)-1) );
		// check if first base pair of seedY is enclosed base pair of seedX
		auto bpX = std::find( xBP.basePairs.begin(), xBP.basePairs.end(), energy.getBasePair(y1,y2) );
		if ( bpX == xBP.basePairs.end())
		{
			return false;
		}
		size_t overlap = 1;
		// collect all base pairs of seedY
		Interaction yBP( *(xBP.s1), *(xBP.s2) ); // use seq info of xBP
		yBP.basePairs.push_back( energy.getBasePair(y1,y2) );
		traceBackSeed( yBP, y1, y2 );
		yBP.basePairs.push_back( energy.getBasePair( y1 + getSeedLength1(y1,y2) - 1, y2 + getSeedLength2(y1,y2)-1) );
		// check if all subsequent bps of  are equal
		auto bpY = yBP.basePairs.begin();
		while( (++bpX) != xBP.basePairs.end()
				// ensure there is still a bp in seedY to compare to
				&& (++bpY) != yBP.basePairs.end()
				// and that the base pairs are equal
				&& *bpX == *bpY )
		{
			overlap++;
		}
		// ensure seedX was fully processed
		// and seeds are overlapping in at least 2 base pairs, ie. the last
		// loop of seedX
		return bpX == xBP.basePairs.end() && overlap > 1;
	}
	// not valid
	return false;
}


//////////////////////////////////////////////////////////////////////////

} // namespace


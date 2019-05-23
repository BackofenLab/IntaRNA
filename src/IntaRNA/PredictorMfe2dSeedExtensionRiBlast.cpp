
#include "IntaRNA/PredictorMfe2dSeedExtensionRiBlast.h"

namespace IntaRNA {

//////////////////////////////////////////////////////////////////////////

PredictorMfe2dSeedExtensionRiBlast::
PredictorMfe2dSeedExtensionRiBlast(
		const InteractionEnergy & energy
		, OutputHandler & output
		, PredictionTracker * predTracker
		, SeedHandler * seedHandlerInstance )
 :
	PredictorMfe2d(energy,output,predTracker)
	, seedHandler(seedHandlerInstance)
	, hybridE_right( 0,0 )
{
	assert( seedHandler.getConstraint().getBasePairs() > 1 );
}

//////////////////////////////////////////////////////////////////////////

PredictorMfe2dSeedExtensionRiBlast::
~PredictorMfe2dSeedExtensionRiBlast()
{
}

//////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dSeedExtensionRiBlast::
predict( const IndexRange & r1, const IndexRange & r2
		, const OutputConstraint & outConstraint )
{
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) <<"predicting mfe interactions with seed in O(n^2) space and O(n^4) time..."; }
	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

	// suboptimal setup check
	if (outConstraint.reportMax>1 && outConstraint.reportOverlap != OutputConstraint::ReportOverlap::OVERLAP_BOTH) {
		throw std::runtime_error("PredictorMfe2dSeedExtensionRiBlast : the enumeration of non-overlapping suboptimal interactions is not supported in this prediction mode");
	}

#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!(r1.isAscending() && r2.isAscending()) )
		throw std::runtime_error("PredictorMfe2d::predict("+toString(r1)+","+toString(r2)+") is not sane");
#endif

	// setup index offset
	energy.setOffset1(r1.from);
	energy.setOffset2(r2.from);
	seedHandler.setOffset1(r1.from);
	seedHandler.setOffset2(r2.from);

	const size_t interaction_size1 = std::min( energy.size1()
			, (r1.to==RnaSequence::lastPos?energy.size1()-1:r1.to)-r1.from+1 );
	const size_t interaction_size2 = std::min( energy.size2()
			, (r2.to==RnaSequence::lastPos?energy.size2()-1:r2.to)-r2.from+1 );

	// compute seed interactions for whole range
	// and check if any seed possible
	if (seedHandler.fillSeed( 0, interaction_size1-1, 0, interaction_size2-1 ) == 0) {
		// trigger empty interaction reporting
		initOptima(outConstraint);
		reportOptima(outConstraint);
		// stop computation
		return;
	}

	// initialize mfe interaction for updates
	initOptima( outConstraint );

	size_t si1 = RnaSequence::lastPos, si2 = RnaSequence::lastPos;
	while( seedHandler.updateToNextSeed(si1,si2
			, 0, interaction_size1+1-seedHandler.getConstraint().getBasePairs()
			, 0, interaction_size2+1-seedHandler.getConstraint().getBasePairs()) )
	{
		E_type seedE = seedHandler.getSeedE(si1, si2);
		size_t sl1 = seedHandler.getSeedLength1(si1, si2);
		size_t sl2 = seedHandler.getSeedLength2(si1, si2);
		const size_t sj1 = si1+sl1-1;
		const size_t sj2 = si2+sl2-1;
		// check if seed fits into interaction range
		if (sj1 > interaction_size1 || sj2 > interaction_size2)
			continue;

		ExtendedSeed extension;
		extension.i1 = si1;
		extension.i2 = si2;
		extension.j1 = sj1;
		extension.j2 = sj2;
		extension.energy = seedE;
		size_t maxMatrixLen1 = energy.getAccessibility1().getMaxLength()-sl1+1;
		size_t maxMatrixLen2 = energy.getAccessibility2().getMaxLength()-sl2+1;

		parallelExtension(extension, std::min(interaction_size1-extension.j1, maxMatrixLen1), std::min(interaction_size2-extension.j2, maxMatrixLen2) );

		// update constraints
		sl1 = extension.j1 - extension.i1;
		sl2 = extension.j2 - extension.i2;
    maxMatrixLen1 = energy.getAccessibility1().getMaxLength()-sl1;
		maxMatrixLen2 = energy.getAccessibility2().getMaxLength()-sl2;

		// EL
		hybridE_left.resize( std::min(extension.i1+1, maxMatrixLen1), std::min(extension.i2+1, maxMatrixLen2) );
		for (size_t i = 0; i < hybridE_left.size1(); i++) {
			for (size_t j = 0; j < hybridE_left.size2(); j++) {
				hybridE_left(i, j) = E_INF;
			}
		}
		fillHybridE_left(extension.i1, extension.i2, outConstraint);

		// ER
		hybridE_right.resize( std::min(interaction_size1-extension.j1, maxMatrixLen1), std::min(interaction_size2-extension.j2, maxMatrixLen2) );
		for (size_t i = 0; i < hybridE_right.size1(); i++) {
			for (size_t j = 0; j < hybridE_right.size2(); j++) {
				hybridE_right(i, j) = E_INF;
			}
		}
		fillHybridE_right(extension.j1, extension.j2, outConstraint);

		// update Optimum for all boundary combinations
		for (int i1 = 0; i1 < hybridE_left.size1(); i1++) {
			// ensure max interaction length in seq 1
			for (int j1 = 0; j1 < hybridE_right.size1() ; j1++) {
				if (extension.j1+j1-extension.i1+i1 > energy.getAccessibility1().getMaxLength()) continue;
				for (int i2 = 0; i2 < hybridE_left.size2(); i2++) {
					if (E_isINF(hybridE_left(i1,i2))) continue;
					// ensure max interaction length in seq 2
					for (int j2 = 0; j2 < hybridE_right.size2() ; j2++) {
						if (extension.j2+j2-extension.i2+i2 > energy.getAccessibility2().getMaxLength()) continue;
						if (E_isINF(hybridE_right(j1,j2))) continue;
						PredictorMfe::updateOptima( extension.i1-i1, extension.j1+j1, extension.i2-i2, extension.j2+j2, extension.energy + hybridE_left(i1,i2) + hybridE_right(j1,j2), true );
					} // j2
				} // i2
			} // j1
		} // i1

	} // si1 / si2

	// report mfe interaction
	reportOptima( outConstraint );

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dSeedExtensionRiBlast::
parallelExtension( PredictorMfe2dSeedExtensionRiBlast::ExtendedSeed & seed
	  	, const size_t max_extension1, size_t max_extension2
	  	)
{

	size_t i1min = seed.i1;
	size_t i2min = seed.i2;

	// extend left
	while (seed.i1 > 1 && seed.i2 > 1) {
		// todo acc1/acc2 maxLength() termination
		if (energy.areComplementary(seed.i1-1,seed.i2-1)) {
			E_type newEnergy = seed.energy + energy.getE_interLeft(seed.i1-1,i1min,seed.i2-1,i2min);
			if (newEnergy < seed.energy) {
				seed.energy = newEnergy;
				i1min = seed.i1-1;
				i2min = seed.i2-1;
			}
		}
		if (i1min-(seed.i1-1) >= parallelDropOutLength) {
			break;
		}
		seed.i1--;
		seed.i2--;
	}
	seed.i1 = i1min;
	seed.i2 = i2min;

	size_t j1min = seed.j1;
	size_t j2min = seed.j2;

	// extend right
	while (seed.j1 < max_extension1-1 && seed.j2 < max_extension2-1) {
		// todo acc1/acc2 maxLength() termination
		if (energy.areComplementary(seed.j1+1,seed.j2+1)) {
			E_type newEnergy = seed.energy + energy.getE_interLeft(j1min,seed.j1+1,j2min,seed.j2+1);
			if (newEnergy < seed.energy) {
				seed.energy = newEnergy;
				j1min = seed.j1+1;
				j2min = seed.j2+1;
			}
		}
		if (seed.j1+1-i1min >= parallelDropOutLength) {
			break;
		}
		seed.j1++;
		seed.j2++;
	}
	seed.j1 = j1min;
	seed.j2 = j2min;

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dSeedExtensionRiBlast::
fillHybridE_left( const size_t j1, const size_t j2
			, const OutputConstraint & outConstraint )
{

	// global vars to avoid reallocation
	size_t i1,i2,k1,k2;
	//////////  FIRST ROUND : COMPUTE HYBRIDIZATION ENERGIES ONLY  ////////////

	// current minimal value
	E_type curMinE = E_INF;
	size_t length = 1;
	while (length < 1+std::min(hybridE_left.size1(), hybridE_left.size2())) {
	  for (i1 = 0; i1 < length; i1++) {
			i2 = length-i1-1;

			// init current cell (e_init if just left (i1,i2) base pair)
			hybridE_left(i1,i2) = i1==0 && i2==0 ? energy.getE_init() : E_INF;

			// check if complementary
			if( i1>0 && i2>0 && energy.areComplementary(j1-i1,j2-i2) ) {
				curMinE = E_INF;

				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				for (k1=i1; k1-- > 0; ) {
					// ensure maximal loop length
					if (i1-k1 > energy.getMaxInternalLoopSize1()+1) break;
				for (k2=i2; k2-- > 0; ) {
					// ensure maximal loop length
					if (i2-k2 > energy.getMaxInternalLoopSize2()+1) break;
					// check if (k1,k2) are valid left boundary
					if ( E_isNotINF( hybridE_left(k1,k2) ) ) {
						curMinE = std::min( curMinE,
								(energy.getE_interLeft(j1-i1,j1-k1,j2-i2,j2-k2)
										+ hybridE_left(k1,k2) )
								);
					}
				} // k2
			  } // k1

				// store value
				hybridE_left(i1,i2) = curMinE;
			}
		}
		if (length >= thoroughDropOutLength) {
			break;
		}
		length++;
	}

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dSeedExtensionRiBlast::
fillHybridE_right( const size_t i1, const size_t i2
			, const OutputConstraint & outConstraint )
{

	// global vars to avoid reallocation
	size_t j1,j2,k1,k2;
	//////////  FIRST ROUND : COMPUTE HYBRIDIZATION ENERGIES ONLY  ////////////

	// current minimal value
	E_type curMinE = E_INF;
	size_t length = 1;
	while (length < 1+std::min(hybridE_right.size1(), hybridE_right.size2())) {
	  for (j1 = 0; j1 < length; j1++) {
			j2 = length-j1-1;

			// init current cell (0 if just left (i1,i2) base pair)
			hybridE_right(j1,j2) = j1==0 && j2==0 ? 0 : E_INF;

			// check if complementary
			if( j1>i1 && j2>i2 && energy.areComplementary(i1+j1,i2+j2) ) {
				curMinE = E_INF;

				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				for (k1=j1; k1-- > 0; ) {
					// ensure maximal loop length
					if (j1-k1 > energy.getMaxInternalLoopSize1()+1) break;
				for (k2=j2; k2-- > 0; ) {
					// ensure maximal loop length
					if (j2-k2 > energy.getMaxInternalLoopSize2()+1) break;
					// check if (k1,k2) are valid left boundary
					if ( E_isNotINF( hybridE_right(k1,k2) ) ) {
						curMinE = std::min( curMinE,
								(energy.getE_interLeft(i1+k1,i1+j1,i2+k2,i2+j2)
										+ hybridE_right(k1,k2) )
								);
					}
				} // k2
			  } // k1

				// store value
				hybridE_right(j1,j2) = curMinE;
			}
		}
		if (length >= thoroughDropOutLength) {
			break;
		}
		length++;
	}

}

////////////////////////////////////////////////////////////////////////////


void
PredictorMfe2dSeedExtensionRiBlast::
traceBack( Interaction & interaction, const OutputConstraint & outConstraint  )
{
	// check if something to trace
	if (interaction.basePairs.size() < 2) {
		return;
	}

#if INTARNA_IN_DEBUG_MODE
	// sanity checks
	if ( interaction.basePairs.size() != 2 ) {
		throw std::runtime_error("PredictorMfe2dSeedExtensionRiBlast::traceBack() : given interaction does not contain boundaries only");
	}
#endif

	// check for single interaction
	if (interaction.basePairs.at(0).first == interaction.basePairs.at(1).first) {
		// delete second boundary (identical to first)
		interaction.basePairs.resize(1);
		// update done
		return;
	}

#if INTARNA_IN_DEBUG_MODE
	// sanity checks
	if ( ! interaction.isValid() ) {
		throw std::runtime_error("PredictorMfe2dSeedExtensionRiBlast::traceBack() : given interaction not valid");
	}
#endif

	// ensure sorting
	interaction.sort();
	// get indices in hybridE for boundary base pairs
	size_t	i1 = energy.getIndex1(interaction.basePairs.at(0)),
			j1 = energy.getIndex1(interaction.basePairs.at(1)),
			i2 = energy.getIndex2(interaction.basePairs.at(0)),
			j2 = energy.getIndex2(interaction.basePairs.at(1))
			;

#if INTARNA_IN_DEBUG_MODE
	// check if intervals are larger enough to contain a seed
	if (std::min(j1-i1,j2-i2)+1 < seedHandler.getConstraint().getBasePairs()) {
		// no seed possible, abort computation
		throw std::runtime_error("PredictorMfe2dSeedExtensionRiBlast::traceBack() : given boundaries "+toString(interaction)+" can not hold a seed of "+toString(seedHandler.getConstraint().getBasePairs())+" base pairs");
	}
#endif

	E_type fullE = interaction.energy;

	size_t si1 = RnaSequence::lastPos, si2 = RnaSequence::lastPos;
	while( seedHandler.updateToNextSeed(si1,si2
			, i1,j1+1-seedHandler.getConstraint().getBasePairs()
			, i2,j2+1-seedHandler.getConstraint().getBasePairs() ) )
	{
		E_type seedE = seedHandler.getSeedE(si1, si2);

		size_t sl1 = seedHandler.getSeedLength1(si1, si2);
		size_t sl2 = seedHandler.getSeedLength2(si1, si2);
		size_t sj1 = si1+sl1-1;
		size_t sj2 = si2+sl2-1;

		ExtendedSeed extension;
		extension.i1 = si1;
		extension.i2 = si2;
		extension.j1 = sj1;
		extension.j2 = sj2;
		extension.energy = seedE;

		parallelExtension(extension, j1+1, j2+1);

		sl1 = extension.j1 - extension.i1;
		sl2 = extension.j2 - extension.i2;
		sj1 = extension.j1;
		sj2 = extension.j2;
		seedE = extension.energy;
		const size_t maxMatrixLen1 = energy.getAccessibility1().getMaxLength()-sl1;
		const size_t maxMatrixLen2 = energy.getAccessibility2().getMaxLength()-sl2;

		// traceback for full extension
		hybridE_left.resize( std::min(extension.i1+1, maxMatrixLen1), std::min(extension.i2+1, maxMatrixLen2) );
		for (size_t i = 0; i < hybridE_left.size1(); i++) {
			for (size_t j = 0; j < hybridE_left.size2(); j++) {
				hybridE_left(i, j) = E_INF;
			}
		}
		fillHybridE_left(extension.i1, extension.i2, outConstraint);
		hybridE_right.resize( std::min(j1-extension.j1+1, maxMatrixLen1), std::min(j2-extension.j2+1, maxMatrixLen2) );
		for (size_t i = 0; i < hybridE_right.size1(); i++) {
			for (size_t j = 0; j < hybridE_right.size2(); j++) {
				hybridE_right(i, j) = E_INF;
			}
		}
		fillHybridE_right(extension.j1, extension.j2, outConstraint);

		if ( E_equal( fullE,
				(energy.getE(i1, j1, i2, j2, seedE + hybridE_left( extension.i1-i1, extension.i2-i2 ) + hybridE_right( j1-sj1, j2-sj2 )))))
		{
			// found seed -> traceback
			// the currently traced value for i1-si1, i2-si2
			E_type curE = hybridE_left(extension.i1-i1,extension.i2-i2);

			// trace back left
			while( i1 != extension.i1 ) {

				// check if just internal loop
				if ( E_equal( curE, (energy.getE_interLeft(i1,extension.i1,i2,extension.i2) + hybridE_left(0,0)) ) )
				{
					break;
				}
				// check all interval splits
				if ( (extension.i1-i1) > 1 && (extension.i2-i2) > 1) {
					// temp variables
					size_t k1,k2;
					bool traceNotFound = true;
					// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
					for (k1=std::min(extension.i1-1,i1+energy.getMaxInternalLoopSize1()+1); traceNotFound && k1>i1; k1--) {
					for (k2=std::min(extension.i2-1,i2+energy.getMaxInternalLoopSize2()+1); traceNotFound && k2>i2; k2--) {
						// check if (k1,k2) are valid left boundary
						if ( E_isNotINF( hybridE_left(extension.i1-k1,extension.i2-k2) ) ) {
							if ( E_equal( curE,
									(energy.getE_interLeft(i1,k1,i2,k2) + hybridE_left(extension.i1-k1,extension.i2-k2)) ) )
							{
						    // stop searching
							  traceNotFound = false;
							  // store splitting base pair
							  interaction.basePairs.push_back( energy.getBasePair(k1,k2) );
							  // trace right part of split
							  i1=k1;
						  	i2=k2;
						  	curE = hybridE_left(extension.i1-i1,extension.i2-i2);
							}
						}
					}
					}
				}

			 } // traceback left

			// trace seed
			size_t i = extension.i1;
			size_t j = extension.i2;
			while (i < extension.j1 && j < extension.j2) {
				if (i > i1 && i < j1 && j > i2 && j < j2) {
				  interaction.basePairs.push_back( energy.getBasePair(i, j) );
				 }
				i++;
				j++;
			}

			// the currently traced value for sj1-j1, sj2-j2
			curE = hybridE_right(j1-sj1,j2-sj2);

			// trace back right
			while( j1 != sj1 ) {

				// check if just internal loop
				if ( E_equal( curE, (energy.getE_interLeft(sj1,j1,sj2,j2) + hybridE_right(0,0)) ) )
				{
					break;
				}
				// check all interval splits
				if ( (j1-sj1) > 1 && (j2-sj2) > 1) {
					// temp variables
					size_t k1,k2;
					bool traceNotFound = true;
					// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
					for (k1=std::min(j1-1,sj1+energy.getMaxInternalLoopSize1()+1); traceNotFound && k1>sj1; k1--) {
					for (k2=std::min(j2-1,sj2+energy.getMaxInternalLoopSize2()+1); traceNotFound && k2>sj2; k2--) {
						// check if (k1,k2) are valid left boundary
						if ( E_isNotINF( hybridE_right(k1-sj1,k2-sj2) ) ) {
							if ( E_equal( curE,
									(energy.getE_interLeft(k1,j1,k2,j2) + hybridE_right(k1-sj1,k2-sj2)) ) )
							{
								// stop searching
								traceNotFound = false;
					   		// store splitting base pair
								interaction.basePairs.push_back( energy.getBasePair(k1,k2) );
								// trace right part of split
								j1=k1;
								j2=k2;
								curE = hybridE_right(j1-sj1,j2-sj2);
							}
						}
					}
					}
				}
			}  // traceback right

			interaction.sort();
			seedHandler.addSeeds( interaction );

			// stop searching for seeds
			return;
		}

  } // si1 / si2

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dSeedExtensionRiBlast::
getNextBest( Interaction & curBest )
{
	INTARNA_NOT_IMPLEMENTED("PredictorMfe2dSeedExtensionRiBlast::getNextBest() not implemented yet");
}

//////////////////////////////////////////////////////////////////////////


} // namespace

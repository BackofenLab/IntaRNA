
#include "IntaRNA/OutputHandlerText.h"

#include <sstream>
#include <iomanip>
#include <numeric>
#include <algorithm>

#if INTARNA_MULITHREADING
	#include <omp.h>
#endif

namespace IntaRNA {

////////////////////////////////////////////////////////////////////////////

OutputHandlerText::
OutputHandlerText(
		const OutputConstraint & outConstraint,
		std::ostream & out,
		const InteractionEnergy & energy,
		const size_t flankingLength_,
		const bool detailedOutput
		)
 :	OutputHandler(outConstraint)
	, out(out)
	, energy(energy)
	, flankingLength(flankingLength_)
	, detailedOutput(detailedOutput)
{
	if (flankingLength < 10) {
		throw std::runtime_error("OutputHandlerText::OutputHandlerText : flankingLength ("+toString(flankingLength_)+") has to be at least 9");
	}
}

////////////////////////////////////////////////////////////////////////////

OutputHandlerText::
~OutputHandlerText()
{

	// special handling if no base pairs present
	if (reportedInteractions == 0) {
		// ensure outputs do not intervene
		std::stringstream outTmp;
		outTmp <<"\n"
			<<"no favorable interaction for "
			<<energy.getAccessibility1().getSequence().getId()
			<<" and "
			<<energy.getAccessibility2().getSequence().getId()
			<<'\n';
#if INTARNA_MULITHREADING
		#pragma omp critical(intarna_omp_outputStreamUpdate)
#endif
		{
			out <<outTmp.str();
		} // omp critical(intarna_omp_outputStreamUpdate)
	}

	out.flush();
}

////////////////////////////////////////////////////////////////////////////

void
OutputHandlerText::
add( const Interaction & i )
{
#if INTARNA_IN_DEBUG_MODE
	// debug checks
	if ( i.basePairs.size() > 0 && ! i.isValid() ) {
		throw std::runtime_error("OutputHandlerText::add() : given interaction is not valid : "+toString(i));
	}
#endif

	// count the interaction
	reportedInteractions++;

	// get interaction start/end per sequence
	const size_t i1 = i.basePairs.begin()->first;
	const size_t j1 = i.basePairs.rbegin()->first;
	const size_t i2 = i.basePairs.begin()->second;
	const size_t j2 = i.basePairs.rbegin()->second;

	const IndexRangeList seedRanges1 = i.getSeedRanges1();

	// decompose printing into several rows
	std::ostringstream s1Unbound;
	std::ostringstream s1Bound;
	std::ostringstream pairing;
	std::ostringstream s2Bound;
	std::ostringstream s2Unbound;

	// left flanking
	// unbound region s1
	s1Unbound.width(flankingLength+3);
	s1Unbound <<std::right;
	if (i1 < flankingLength) {
		// full sequence prefix
		s1Unbound <<("5'-" + i.s1->asString().substr( 0, i1));
	} else {
		// prefix shortening needed
		s1Unbound <<("5'-"
					+ i.s1->asString().substr( 0,3)
					+ "..."
					+ i.s1->asString().substr( (size_t)std::max(0,(int)i1+6-(int)flankingLength), flankingLength-6)
					);
	}
	// bound region
	s1Bound.width(flankingLength+3); s1Bound <<' ';
	pairing.width(flankingLength+3); pairing <<' ';
	s2Bound.width(flankingLength+3); s2Bound <<' ';
	// unbound region s2
	s2Unbound.width(flankingLength+3);
	s2Unbound <<std::right; // <<("3'-" + i.s2->asString().substr( (size_t)std::max(0,(int)i2-(int)flankingLength), std::min(i2,flankingLength) ));
//	if (i2 < flankingLength) {
//		// full sequence prefix
//		s2Unbound <<("3'-" + i.s2->asString().substr( 0, i2));
//	} else {
//		// prefix shortening needed
//		s2Unbound <<("3'-"
//					+ i.s2->asString().substr( 0,3)
//					+ "..."
//					+ i.s2->asString().substr( (size_t)std::max(0,(int)i2+6-(int)flankingLength), flankingLength-6)
//					);
//	}
	if (i2+flankingLength > i.s2->size()) {
		// add remaining sequence
		s2Unbound <<("3'-" + reverse(i.s2->asString().substr(i2+1)));
	} else {
		// shortening needed
		s2Unbound <<("3'-"
					+ reverse(i.s2->asString().substr(i.s2->size()-3))
					+ "..."
					+ reverse(i.s2->asString().substr(i2+1, flankingLength-6))
					);
	}

	// interaction start left
	Interaction::PairingVec::const_iterator leftBP = i.basePairs.begin();
	// print start
	s1Unbound.width(1);	 s1Unbound <<std::left <<' ';
	s1Bound.width(1);	 s1Bound   <<std::left <<i.s1->asString().at(leftBP->first);
	pairing.width(1);	 pairing   <<std::left;
	if (seedRanges1.covers(leftBP->first)) {
		pairing   <<'+';
	} else
	if (energy.isGU(energy.getIndex1(*leftBP), energy.getIndex2(*leftBP))) {
		pairing   <<':';
	} else {
		pairing   <<'|';
	}
	s2Bound.width(1);	 s2Bound   <<std::left <<i.s2->asString().at(leftBP->second);
	s2Unbound.width(1);	 s2Unbound <<std::left <<' ';

	// iterate loops in interaction region
	Interaction::PairingVec::const_iterator curBP = i.basePairs.begin();
	size_t loop1=0, loop2=0, loop=0, interactionLength = 1;
	for (++curBP; curBP != i.basePairs.end(); ++curBP, ++leftBP) {
		// handle duplicated BPs (might happen due to explicit seeds containing only a single bp)
		if (*curBP == *leftBP) {continue;}
		// handle internal loop region
		// get specific loop lengths
		loop1 = curBP->first - leftBP->first -1;
		loop2 = leftBP->second - curBP->second -1;
		loop = std::max(loop1,loop2);
		// print unbound loop regions
		if (loop>0) {
			// unbound region s1
			if (loop1 > 0) {
				s1Unbound <<i.s1->asString().substr( leftBP->first +1, loop1 );
				// fill missing positions
				if (loop1 < loop) {
					s1Unbound.width(loop-loop1);
					s1Unbound <<std::setfill('-') <<'-' <<std::setfill(' ');
				}
			} else {
				s1Unbound.width(loop);
				s1Unbound <<std::setfill('-') <<'-' <<std::setfill(' ');
			}
			// bound region
			s1Bound.width(loop); s1Bound <<' ';
			pairing.width(loop); pairing <<' ';
			s2Bound.width(loop); s2Bound <<' ';
			// unbound region s2
			if (loop2 > 0) {
				s2Unbound <<reverse(i.s2->asString().substr( curBP->second +1, loop2 ));
				// fill missing positions
				if (loop2 < loop) {
					s2Unbound.width(loop-loop2);
					s2Unbound <<std::setfill('-') <<'-' <<std::setfill(' ');
				}
			} else {
				s2Unbound.width(loop);
				s2Unbound <<std::setfill('-') <<'-' <<std::setfill(' ');
			}
		}
		interactionLength += loop;

		// print current base pair (right end of internal loop)
		s1Unbound.width(1);	 s1Unbound <<' ';
		s1Bound.width(1);	 s1Bound   <<i.s1->asString().at(curBP->first);
		pairing.width(1);
		if (seedRanges1.covers(curBP->first)) {
			pairing   <<'+';
		} else
		if (energy.isGU(energy.getIndex1(*curBP), energy.getIndex2(*curBP))) {
			pairing   <<':';
		} else {
			pairing   <<'|';
		}
		s2Bound.width(1);	 s2Bound   <<i.s2->asString().at(curBP->second);
		s2Unbound.width(1);	 s2Unbound <<' ';
		interactionLength++;
	}

	// flanking right
	// unbound region s1
	if (j1+flankingLength > i.s1->size()) {
		// add remaining sequence
		s1Unbound <<i.s1->asString().substr(j1+1);
	} else {
		// shortening needed
		s1Unbound <<i.s1->asString().substr(j1+1, flankingLength-6)
					<<"..."
					<<i.s1->asString().substr(i.s1->size()-3);
	}
	s1Unbound <<"-3'";
	// unbound region s2
//	if (j2+flankingLength > i.s2->size()) {
//		// add remaining sequence
//		s2Unbound <<i.s2->asString().substr(j2+1);
//	} else {
//		// shortening needed
//		s2Unbound <<i.s2->asString().substr(j2+1, flankingLength-6)
//					<<"..."
//					<<i.s2->asString().substr(i.s2->size()-3);
//	}
	if (j2 < flankingLength) {
		// full sequence prefix
		s2Unbound <<reverse(i.s2->asString().substr( 0, j2));
	} else {
		// prefix shortening needed
		s2Unbound 	<<reverse(i.s2->asString().substr( (size_t)std::max(0,(int)j2+6-(int)flankingLength), flankingLength-6))
					<<"..."
					<<reverse(i.s2->asString().substr( 0,3))
					;
	}
	s2Unbound <<"-5'";

	// position information (first and last interacting positions in sequence)
	// start indexing with 1 instead of 0
	std::ostringstream pos1, pos1tag;
	std::string i1Str = toString(i.s1->getInOutIndex(i1));
	std::string j1Str = toString(i.s1->getInOutIndex(j1));
	if (i.s1->getInOutIndex(0)<0) {
		// explicit "+" annotation of positive indices in relative indexing
		if (i.s1->getInOutIndex(i1)>0) {
			i1Str = "+" + i1Str;
		}
		if (i.s1->getInOutIndex(j1)>0) {
			j1Str = "+" + j1Str;
		}
	}
	pos1	<<std::setw(flankingLength+3+1) <<std::right <<i1Str;
	pos1tag	<<std::setw(flankingLength+3+1) <<'|';
	// put right (3') end only if not overlapping with 5' start position (left)
	if (i1+1 < j1 && i1Str.size()+2 < flankingLength+3+interactionLength) {
		pos1 <<std::setw(interactionLength -2) <<' '
				<<std::setw(1) <<std::left <<j1Str;
		pos1tag	<<std::setw(interactionLength - 1) <<'|';
	}
	std::ostringstream pos2, pos2tag;
	std::string i2Str = toString(i.s2->getInOutIndex(i2));
	std::string j2Str = toString(i.s2->getInOutIndex(j2));
	if (i.s2->getInOutIndex(0)<0) {
		// explicit "+" annotation of positive indices in relative indexing
		if (i.s2->getInOutIndex(i2)>0) {
			i2Str = "+" + i2Str;
		}
		if (i.s2->getInOutIndex(j2)>0) {
			j2Str = "+" + j2Str;
		}
	}
	// put left (3') end only if not overlapping with 5' start position (right)
	pos2	<<std::setw(flankingLength+3+1);
	pos2tag	<<std::setw(flankingLength+3+1);
	if (i2==j2 || i2Str.size()+2 < (flankingLength+3+interactionLength)) {
		pos2	 <<std::right <<i2Str;
		pos2tag	 <<'|';
	} else {
		pos2	 <<' ';
		pos2tag	 <<' ';
	}
	if (i2 > j2) {
		pos2 <<std::setw(interactionLength - 2) <<' '
				<<std::setw(1) <<std::left <<j2Str;
		pos2tag	<<std::setw(interactionLength - 1) <<'|';
	}

	// get individual energy contributions
	InteractionEnergy::EnergyContributions contr = energy.getE_contributions(i);

	{
		// ensure outputs do not intervene
		std::stringstream outTmp;

		// print full interaction to output stream
		outTmp <<'\n'
			// get ID of s1
			<<i.s1->getId() <<'\n'
			// get position in s1
			<<pos1.str() <<'\n'
			<<pos1tag.str() <<'\n'
			// print collected interaction stuff
			<<s1Unbound.str() <<'\n'
			<<s1Bound.str() <<'\n'
			<<pairing.str() <<'\n'
			<<s2Bound.str() <<'\n'
			<<s2Unbound.str() <<'\n'
			// get position in reversed index order for s2
			<<pos2tag.str() <<'\n'
			<<pos2.str() <<'\n'
			// get ID of s2
			<<i.s2->getId() <<'\n'
			;

		if (detailedOutput) {
			outTmp
				// interaction range
				<<"\n"
				<<"interaction seq1   = "<<i.s1->getInOutIndex(i.basePairs.begin()->first)<<".."<<i.s1->getInOutIndex(i.basePairs.rbegin()->first) <<'\n'
				<<"interaction seq2   = "<<i.s2->getInOutIndex(i.basePairs.rbegin()->second)<<".."<<i.s2->getInOutIndex(i.basePairs.begin()->second) <<'\n'
				;
		} // detailed
			// print energy
		outTmp
			<<"\n"
			<<"interaction energy = "<<E_2_Ekcal(i.energy) <<" kcal/mol\n"
			;

		if (detailedOutput) {
			outTmp
				<<"  = E(init)        = "<<E_2_Ekcal(contr.init)<<'\n'
				<<"  + E(loops)       = "<<E_2_Ekcal(contr.loops)<<'\n'
				<<"  + E(dangleLeft)  = "<<E_2_Ekcal(contr.dangleLeft)<<'\n'
				<<"  + E(dangleRight) = "<<E_2_Ekcal(contr.dangleRight)<<'\n'
				<<"  + E(endLeft)     = "<<E_2_Ekcal(contr.endLeft)<<'\n'
				<<"  + E(endRight)    = "<<E_2_Ekcal(contr.endRight)<<'\n'
				<<"    : E(hybrid)    = "<<E_2_Ekcal((i.energy-contr.ED1-contr.ED2))<<'\n'
				<<"  + ED(seq1)       = "<<E_2_Ekcal(contr.ED1)<<'\n'
				<<"    : Pu(seq1)     = "<<(E_equal(contr.ED2,0) ? Z_type(1) : energy.getBoltzmannWeight(contr.ED1))<<'\n'
				<<"  + ED(seq2)       = "<<E_2_Ekcal(contr.ED2)<<'\n'
				<<"    : Pu(seq2)     = "<<(E_equal(contr.ED2,0) ? Z_type(1) : energy.getBoltzmannWeight(contr.ED2))<<'\n'
				;
			if (!E_equal(contr.energyAdd,E_type(0))) {
				outTmp
				<<"  + E(add)         = "<<E_2_Ekcal(contr.energyAdd)<<'\n'
				;
			}

			// print seed information if available
			if (i.seed != NULL) {
				const std::string listSep = " | ";
				const auto bestSeed = i.seed->begin();
				// print via std::for_each instead of std::accumulate due to rounding issues of boost::lexical_cast or std::to_string
				// since sometimes (float(int)/100.0) gives strings with 10-5 deviations of expected value
				outTmp
					<<"\nseed seq1   = "<<i.s1->getInOutIndex(bestSeed->bp_i.first)<<".."<<i.s1->getInOutIndex(bestSeed->bp_j.first);
				if (!outConstraint.bestSeedOnly)
					std::for_each( ++(i.seed->begin()), i.seed->end(), [&]( const Interaction::Seed & s) {
								 outTmp << listSep <<i.s1->getInOutIndex(s.bp_i.first)<<".."<<i.s1->getInOutIndex(s.bp_j.first);
								});
				outTmp
					<<"\nseed seq2   = "<<i.s2->getInOutIndex(bestSeed->bp_j.second)<<".."<<i.s2->getInOutIndex(bestSeed->bp_i.second);
				if (!outConstraint.bestSeedOnly)
					std::for_each( ++(i.seed->begin()), i.seed->end(), [&]( const Interaction::Seed & s) {
								 outTmp << listSep <<i.s2->getInOutIndex(s.bp_j.second)<<".."<<i.s2->getInOutIndex(s.bp_i.second);
								});
				outTmp
					<<"\nseed energy = "<<E_2_Ekcal(bestSeed->energy);
				if (!outConstraint.bestSeedOnly)
					std::for_each( ++(i.seed->begin()), i.seed->end(), [&]( const Interaction::Seed & s) {
								 outTmp << listSep << E_2_Ekcal(s.energy);
								});
				outTmp
//						<<" kcal/mol"
					<<"\nseed ED1    = "<<E_2_Ekcal(energy.getED1( bestSeed->bp_i.first, bestSeed->bp_j.first ));
				if (!outConstraint.bestSeedOnly)
					std::for_each( ++(i.seed->begin()), i.seed->end(), [&]( const Interaction::Seed & s) {
								 outTmp << listSep << E_2_Ekcal(energy.getED1( s.bp_i.first, s.bp_j.first ));
								});
				outTmp
//						<<" kcal/mol"
					<<"\nseed ED2    = "<<E_2_Ekcal(energy.getAccessibility2().getAccessibilityOrigin().getED( bestSeed->bp_j.second, bestSeed->bp_i.second ));
				if (!outConstraint.bestSeedOnly)
					std::for_each( ++(i.seed->begin()), i.seed->end(), [&]( const Interaction::Seed & s) {
								 outTmp << listSep << E_2_Ekcal(energy.getAccessibility2().getAccessibilityOrigin().getED( s.bp_j.second, s.bp_i.second ));
								});
				outTmp
//						<<" kcal/mol"
					<<"\nseed Pu1    = "<<(E_equal(energy.getED1( bestSeed->bp_i.first, bestSeed->bp_j.first ),0) ? 1 : energy.getBoltzmannWeight(energy.getED1( bestSeed->bp_i.first, bestSeed->bp_j.first )));
				if (!outConstraint.bestSeedOnly)
					std::for_each( ++(i.seed->begin()), i.seed->end(), [&]( const Interaction::Seed & s) {
								 outTmp << listSep <<(E_equal(energy.getED1( s.bp_i.first, s.bp_j.first ),0) ? 1 : energy.getBoltzmannWeight(energy.getED1( s.bp_i.first, s.bp_j.first )));
								});
				outTmp
					<<"\nseed Pu2    = " <<(E_equal(energy.getAccessibility2().getAccessibilityOrigin().getED( bestSeed->bp_j.second, bestSeed->bp_i.second ),0) ? 1 : energy.getBoltzmannWeight(energy.getAccessibility2().getAccessibilityOrigin().getED( bestSeed->bp_j.second, bestSeed->bp_i.second )));
				if (!outConstraint.bestSeedOnly)
					std::for_each( ++(i.seed->begin()), i.seed->end(), [&]( const Interaction::Seed & s) {
								 outTmp << listSep <<(E_equal(energy.getAccessibility2().getAccessibilityOrigin().getED( s.bp_j.second, s.bp_i.second ),0) ? 1 : energy.getBoltzmannWeight(energy.getAccessibility2().getAccessibilityOrigin().getED( s.bp_j.second, s.bp_i.second )));
								});
				outTmp
					<<'\n'
					;
			} // seed
		} // detailed

		// ensemble output if available
		if (!Z_equal(Z,Z_type(0)) && Z_isNotINF(Z)) {
		outTmp
			<<"\n"
			<<"ensemble energy    = "<<E_2_Ekcal(energy.getE(Z)) <<" kcal/mol\n"
			;
		}

		// ensure outputs do not intervene
	#if INTARNA_MULITHREADING
		#pragma omp critical(intarna_omp_outputStreamUpdate)
	#endif
		{
			out <<outTmp.str();
		} // omp critical(intarna_omp_outputStreamUpdate)
	}

}


////////////////////////////////////////////////////////////////////////////

} // namespace

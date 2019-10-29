
#include "IntaRNA/OutputHandlerCsv.h"

#if INTARNA_MULITHREADING
	#include <omp.h>
#endif

#include <boost/algorithm/string.hpp>

namespace IntaRNA {

////////////////////////////////////////////////////////////////////////

std::map<OutputHandlerCsv::ColType,std::string> OutputHandlerCsv::colType2string;

////////////////////////////////////////////////////////////////////////

const std::string OutputHandlerCsv::notAvailable = "NAN";

////////////////////////////////////////////////////////////////////////

const OutputHandlerCsv::ColTypeList OutputHandlerCsv::colTypeNumericSort(
		OutputHandlerCsv::string2list(
		"start1,end1,start2,end2"
		",E,ED1,ED2,Pu1,Pu2,E_init,E_loops,E_dangleL,E_dangleR,E_endL,E_endR,E_hybrid,E_norm,E_hybridNorm,E_add"
		",w"
		",seedStart1,seedEnd1,seedStart2,seedEnd2,seedE,seedED1,seedED2,seedPu1,seedPu2"
		",Eall,Eall1,Eall2"
		",Zall,Zall1,Zall2"
		",Etotal,EallTotal"
		",P_E"
		));

////////////////////////////////////////////////////////////////////////

OutputHandlerCsv::OutputHandlerCsv(
		const OutputConstraint & outConstraint
		,  std::ostream & out_
		, const InteractionEnergy & energy
		, const ColTypeList columns
		, const std::string& colSep
		, const bool printHeader
		, const std::string& listSep
		)
 :	OutputHandler(outConstraint)
	, out(out_)
	, energy(energy)
	, columns(columns)
	, colSep(colSep)
	, listSep(listSep)
{
	// init mapping of coltypes to string
	initColType2string();

	// print CSV header of column names
	if (printHeader) {
		// ensure outputs do not intervene
#if INTARNA_MULITHREADING
		#pragma omp critical(intarna_omp_outputStreamUpdate)
#endif
		{
			out <<getHeader(columns,colSep);
		}
	}
}

////////////////////////////////////////////////////////////////////////

OutputHandlerCsv::~OutputHandlerCsv()
{
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_outputStreamUpdate)
#endif
	{
		// force output
		out.flush();
	}
}

////////////////////////////////////////////////////////////////////////

void
OutputHandlerCsv::
add( const Interaction & i )
{
#if INTARNA_IN_DEBUG_MODE
	// debug checks
	if ( i.basePairs.size() > 0 && ! i.isValid() ) {
		throw std::runtime_error("OutputHandlerCsv::add() : given interaction is not valid : "+toString(i));
	}
#endif

	// special handling if no base pairs present
	if (i.basePairs.size() == 0) {
		return;
	}

	// count the interaction
	reportedInteractions++;

	// get interaction start/end per sequence
	const size_t i1 = i.basePairs.begin()->first;
	const size_t j1 = i.basePairs.rbegin()->first;
	const size_t i2 = i.basePairs.begin()->second;
	const size_t j2 = i.basePairs.rbegin()->second;

	// get individual energy contributions
	InteractionEnergy::EnergyContributions contr = energy.getE_contributions(i);

	// ensure outputs do not intervene
	{
		std::stringstream outTmp;

		for (auto col = columns.begin(); col != columns.end(); col++) {
			// print separator if needed
			if (col != columns.begin()) {
				outTmp <<colSep;
			}
			// print this column information
			switch ( *col ) {

			case id1:
				// ensure no colSeps are contained
				outTmp <<boost::replace_all_copy(i.s1->getId(), colSep, "_");
				break;

			case id2:
				// ensure no colSeps are contained
				outTmp <<boost::replace_all_copy(i.s2->getId(), colSep, "_");
				break;

			case seq1:
				outTmp <<i.s1->asString();
				break;

			case seq2:
				outTmp <<i.s2->asString();
				break;

			case subseq1:
				outTmp <<i.s1->asString().substr(i1, j1-i1+1);
				break;

			case subseq2:
				outTmp <<i.s2->asString().substr(j2, i2-j2+1);
				break;

			case subseqDP:
				outTmp <<i.s1->asString().substr(i1, j1-i1+1)
					<<'&'
					<<i.s2->asString().substr(j2, i2-j2+1);
				break;

			case subseqDB:
				outTmp <<i.s1->getInOutIndex(i1)
					<<i.s1->asString().substr(i1, j1-i1+1)
					<<'&'
					<<i.s2->getInOutIndex(j2)
					<<i.s2->asString().substr(j2, i2-j2+1);
				break;

			case start1:
				outTmp <<i.s1->getInOutIndex(i1);
				break;

			case end1:
				outTmp <<i.s1->getInOutIndex(j1);
				break;

			case start2:
				outTmp <<i.s2->getInOutIndex(j2);
				break;

			case end2:
				outTmp <<i.s2->getInOutIndex(i2);
				break;

			case hybridDP:
				outTmp <<Interaction::dotBracket( i );
				break;

			case hybridDB:
				outTmp <<Interaction::dotBar( i );
				break;

			case hybridDPfull:
				outTmp <<Interaction::dotBracket( i, '(', ')', true );
				break;

			case hybridDBfull:
				outTmp <<Interaction::dotBar( i, true );
				break;

			case bpList: {
				auto bp = i.basePairs.begin();
				outTmp <<'(' <<i.s1->getInOutIndex(bp->first) <<',' <<i.s2->getInOutIndex(bp->second) <<')' ;
				for( ++bp; bp != i.basePairs.end(); bp++ ) {
					outTmp << listSep <<'(' <<i.s1->getInOutIndex(bp->first) <<',' <<i.s2->getInOutIndex(bp->second) <<')';
				};
				break;}

			case E:
				outTmp <<E_2_Ekcal(i.energy);
				break;

			case Etotal:
				if ( E_isINF(energy.getEall1()) || E_isINF(energy.getEall2()) )
					outTmp << notAvailable;
				else
					outTmp <<E_2_Ekcal( i.energy + energy.getEall1() + energy.getEall2() );
				break;

			case ED1:
				outTmp <<E_2_Ekcal(contr.ED1);
				break;

			case ED2:
				outTmp <<E_2_Ekcal(contr.ED2);
				break;

			case Pu1:
				outTmp <<(E_equal(contr.ED1,0) ? Z_type(1) : energy.getBoltzmannWeight(contr.ED1));
				break;

			case Pu2:
				outTmp <<(E_equal(contr.ED2,0) ? Z_type(1) : energy.getBoltzmannWeight(contr.ED2));
				break;

			case E_init:
				outTmp <<E_2_Ekcal(contr.init);
				break;

			case E_loops:
				outTmp <<E_2_Ekcal(contr.loops);
				break;

			case E_dangleL:
				outTmp <<E_2_Ekcal(contr.dangleLeft);
				break;

			case E_dangleR:
				outTmp <<E_2_Ekcal(contr.dangleRight);
				break;

			case E_endL:
				outTmp <<E_2_Ekcal(contr.endLeft);
				break;

			case E_endR:
				outTmp <<E_2_Ekcal(contr.endRight);
				break;

			case E_hybrid:
				outTmp <<E_2_Ekcal(i.energy - contr.ED1 - contr.ED2);
				break;

			case E_norm:
				outTmp <<E_2_Ekcal(i.energy) / std::log( energy.size1() * energy.size2() );
				break;

			case E_hybridNorm:
				outTmp <<E_2_Ekcal(i.energy - contr.ED1 - contr.ED2) / std::log( energy.size1() * energy.size2() );
				break;

			case E_add:
				outTmp <<E_2_Ekcal(contr.energyAdd);
				break;

			case w:
				outTmp <<energy.getBoltzmannWeight(i.energy);
				break;

			case seedStart1:
				if (i.seed == NULL) {
					outTmp <<notAvailable;
				} else {
					if (!outConstraint.bestSeedOnly) {
						// generate list
						outTmp << std::accumulate( std::next( i.seed->begin() )
										, i.seed->end()
										, toString( i.s1->getInOutIndex(i.seed->begin()->bp_i.first) ) // start with first element
										, [&](std::string a, Interaction::Seed s) {
											 return std::move(a) + listSep + toString( i.s1->getInOutIndex(s.bp_i.first) ); // extend list
										 });
					} else {
						outTmp <<i.s1->getInOutIndex(i.seed->begin()->bp_i.first);
					}
				}
				break;

			case seedEnd1:
				if (i.seed == NULL) {
					outTmp <<notAvailable;
				} else {
					if (!outConstraint.bestSeedOnly) {
						// generate list
						outTmp << std::accumulate( std::next( i.seed->begin() )
										, i.seed->end()
										, toString( i.s1->getInOutIndex(i.seed->begin()->bp_j.first) ) // start with first element
										, [&](std::string a, Interaction::Seed s) {
											 return std::move(a) + listSep + toString( i.s1->getInOutIndex(s.bp_j.first) ); // extend list
										 });
					} else {
						outTmp <<i.s1->getInOutIndex(i.seed->begin()->bp_j.first);
					}
				}
				break;

			case seedStart2:
				if (i.seed == NULL) {
					outTmp <<notAvailable;
				} else {
					if (!outConstraint.bestSeedOnly) {
						// generate list
						outTmp << std::accumulate( std::next( i.seed->begin() )
										, i.seed->end()
										, toString( i.s2->getInOutIndex(i.seed->begin()->bp_j.second) ) // start with first element
										, [&](std::string a, Interaction::Seed s) {
											 return std::move(a) + listSep + toString( i.s2->getInOutIndex(s.bp_j.second) ); // extend list
										 });
					} else {
						outTmp <<i.s2->getInOutIndex(i.seed->begin()->bp_j.second);
					}
				}
				break;

			case seedEnd2:
				if (i.seed == NULL) {
					outTmp <<notAvailable;
				} else {
					if (!outConstraint.bestSeedOnly) {
						// generate list
						outTmp << std::accumulate( std::next( i.seed->begin() )
										, i.seed->end()
										, toString( i.s2->getInOutIndex(i.seed->begin()->bp_i.second) ) // start with first element
										, [&](std::string a, Interaction::Seed s) {
											 return std::move(a) + listSep + toString( i.s2->getInOutIndex(s.bp_i.second) ); // extend list
										 });
					} else {
						outTmp <<i.s2->getInOutIndex(i.seed->begin()->bp_i.second);
					}
				}
				break;

			case seedE:
				if (i.seed == NULL) {
					outTmp <<notAvailable;
				} else {
					outTmp << E_2_Ekcal(i.seed->begin()->energy);
					if (!outConstraint.bestSeedOnly) {
						// generate list
						// print via std::for_each instead of std::accumulate due to rounding issues of boost::lexical_cast or std::to_string
						// since sometimes (float(int)/100.0) gives strings with 10E-5 deviations of expected value
						std::for_each( ++(i.seed->begin()), i.seed->end(), [&]( const Interaction::Seed & s) {
									 outTmp << listSep << E_2_Ekcal(s.energy);
									});
					}
				}
				break;

			case seedED1:
				if (i.seed == NULL) {
					outTmp <<notAvailable;
				} else {
					outTmp << E_2_Ekcal(energy.getED1( i.seed->begin()->bp_i.first, i.seed->begin()->bp_j.first ));
					if (!outConstraint.bestSeedOnly) {
						// generate list
						std::for_each( ++(i.seed->begin()), i.seed->end(), [&]( const Interaction::Seed & s) {
									 outTmp << listSep << E_2_Ekcal(energy.getED1( s.bp_i.first, s.bp_j.first ));
									});
					}
				}
				break;

			case seedED2:
				if (i.seed == NULL) {
					outTmp <<notAvailable;
				} else {
					outTmp << E_2_Ekcal(energy.getAccessibility2().getAccessibilityOrigin().getED( i.seed->begin()->bp_j.second, i.seed->begin()->bp_i.second ));
					if (!outConstraint.bestSeedOnly) {
						// generate list
						std::for_each( ++(i.seed->begin()), i.seed->end(), [&]( const Interaction::Seed & s) {
							outTmp << listSep << E_2_Ekcal(energy.getAccessibility2().getAccessibilityOrigin().getED( s.bp_j.second, s.bp_i.second ));
						});
					}
				}
				break;

			case seedPu1:
				if (i.seed == NULL) {
					outTmp <<notAvailable;
				} else {
					if (!outConstraint.bestSeedOnly) {
						// generate list
						outTmp << std::accumulate( std::next( i.seed->begin() )
										, i.seed->end()
										, toString( energy.getBoltzmannWeight( energy.getED1( i.seed->begin()->bp_i.first, i.seed->begin()->bp_j.first )) ) // start with first element
										, [&](std::string a, Interaction::Seed s) {
											 return std::move(a) + listSep + toString( energy.getBoltzmannWeight( energy.getED1( s.bp_i.first, s.bp_j.first ) ) ); // extend list
										 });
					} else {
						outTmp <<energy.getBoltzmannWeight( energy.getED1( i.seed->begin()->bp_i.first, i.seed->begin()->bp_j.first ));
					}
				}
				break;

			case seedPu2:
				if (i.seed == NULL) {
					outTmp <<notAvailable;
				} else {
					if (!outConstraint.bestSeedOnly) {
						// generate list
						outTmp << std::accumulate( std::next( i.seed->begin() )
										, i.seed->end()
										, toString( energy.getBoltzmannWeight( energy.getAccessibility2().getAccessibilityOrigin().getED( i.seed->begin()->bp_j.second, i.seed->begin()->bp_i.second )) ) // start with first element
										, [&](std::string a, Interaction::Seed s) {
											 return std::move(a) + listSep + toString( energy.getBoltzmannWeight( energy.getAccessibility2().getAccessibilityOrigin().getED( s.bp_j.second, s.bp_i.second ) ) ); // extend list
										 });
					} else {
						outTmp <<energy.getBoltzmannWeight( energy.getAccessibility2().getAccessibilityOrigin().getED( i.seed->begin()->bp_j.second, i.seed->begin()->bp_i.second ));
					}
				}
				break;

			case Eall:
				if ( Z_equal(Z,Z_type(0)) ) outTmp << notAvailable; else outTmp <<E_2_Ekcal(energy.getE(Z));
				break;

			case Eall1:
				if ( E_isINF(energy.getEall1()) ) outTmp << notAvailable; else outTmp <<E_2_Ekcal(energy.getEall1());
				break;

			case Eall2:
				if ( E_isINF(energy.getEall2()) ) outTmp << notAvailable; else outTmp <<E_2_Ekcal(energy.getEall2());
				break;

			case EallTotal:
				if ( Z_equal(Z,Z_type(0)) || E_isINF(energy.getEall1()) || E_isINF(energy.getEall2()) )
					outTmp << notAvailable;
				else
					outTmp <<E_2_Ekcal( energy.getE(Z) + energy.getEall1() + energy.getEall2() );
				break;

			case Zall:
				if ( Z_equal(Z,Z_type(0)) ) outTmp << notAvailable; else outTmp <<Z;
				break;

			case Zall1:
				if ( E_isINF(energy.getEall1()) ) outTmp << notAvailable; else outTmp <<energy.getBoltzmannWeight(energy.getEall1());
				break;

			case Zall2:
				if ( E_isINF(energy.getEall2()) ) outTmp << notAvailable; else outTmp <<energy.getBoltzmannWeight(energy.getEall2());
				break;

			case P_E:
				if ( Z_equal(Z,Z_type(0)) ) outTmp << notAvailable; else outTmp <<(energy.getBoltzmannWeight(i.energy)/Z);
				break;

			default : throw std::runtime_error("OutputHandlerCsv::add() : unhandled ColType '"+colType2string[*col]+"'");
			}
		}
		outTmp <<'\n';
	#if INTARNA_MULITHREADING
		#pragma omp critical(intarna_omp_outputStreamUpdate)
	#endif
		{
			out << outTmp.str();
		} // omp critical(intarna_omp_outputStreamUpdate)
	}

}

////////////////////////////////////////////////////////////////////////

} // namespace

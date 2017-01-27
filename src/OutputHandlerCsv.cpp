
#include "OutputHandlerCsv.h"

////////////////////////////////////////////////////////////////////////

std::map<OutputHandlerCsv::ColType,std::string> OutputHandlerCsv::colType2string;

////////////////////////////////////////////////////////////////////////

OutputHandlerCsv::OutputHandlerCsv(
		  std::ostream & out
		, const InteractionEnergy & energy
		, const ColTypeList colOrder
		, const std::string& colSep
		, const bool printHeader
		)
 :	out(out)
	, energy(energy)
	, colOrder(colOrder)
	, colSep(colSep)
{
	// init mapping of coltypes to string
	initColType2string();

	// print CSV header of column names
	if (printHeader) {
		out <<getHeader(colOrder,colSep);
	}
}

////////////////////////////////////////////////////////////////////////

OutputHandlerCsv::~OutputHandlerCsv()
{
}

////////////////////////////////////////////////////////////////////////

void
OutputHandlerCsv::
add( const Interaction & i )
{
#if IN_DEBUG_MODE
	// debug checks
	if ( i.basePairs.size() > 0 && ! i.isValid() ) {
		throw std::runtime_error("OutputHandlerCsv::add() : given interaction is not valid : "+toString(i));
	}
#endif

	// special handling if no base pairs present
	if (i.basePairs.size() == 0) {
		return;
	}

	// get interaction start/end per sequence
	const size_t i1 = i.basePairs.begin()->first;
	const size_t j1 = i.basePairs.rbegin()->first;
	const size_t i2 = i.basePairs.begin()->second;
	const size_t j2 = i.basePairs.rbegin()->second;

	// get individual energy contributions
	InteractionEnergy::EnergyContributions contr = energy.getE_contributions(i);

	for (auto col = colOrder.begin(); col != colOrder.end(); col++) {
		// print separator if needed
		if (col != colOrder.begin()) {
			out <<colSep;
		}
		// print this column information
		switch ( *col ) {

		case id1:
			out <<energy.getAccessibility1().getSequence().getId();
			break;

		case id2:
			out <<energy.getAccessibility2().getSequence().getId();
			break;

		case seq1:
			out <<energy.getAccessibility1().getSequence().asString();
			break;

		case seq2:
			out <<energy.getAccessibility2().getAccessibilityOrigin().getSequence().asString();
			break;

		case subseq1:
			out <<energy.getAccessibility1().getSequence().asString().substr(i1, j1-i1+1);
			break;

		case subseq2:
			out <<energy.getAccessibility2().getAccessibilityOrigin().getSequence().asString().substr(j2, i2-j2+1);
			break;

		case subseqDP:
			out <<energy.getAccessibility1().getSequence().asString().substr(i1, j1-i1+1)
				<<'&'
				<<energy.getAccessibility2().getAccessibilityOrigin().getSequence().asString().substr(j2, i2-j2+1);
			break;

		case subseqDB:
			out <<(i1+1)
				<<energy.getAccessibility1().getSequence().asString().substr(i1, j1-i1+1)
				<<'&'
				<<(j2+1)
				<<energy.getAccessibility2().getAccessibilityOrigin().getSequence().asString().substr(j2, i2-j2+1);
			break;

		case hybridDP:
			out <<Interaction::dotBracket( i );
			break;

		case hybridDB:
			out <<Interaction::dotBar( i );
			break;

		case E:
			out <<i.energy;
			break;

		case E_seed:
			if (i.seedRange == NULL) {
				out <<std::numeric_limits<E_type>::signaling_NaN();
			} else {
				out <<i.seedRange->energy;
			}
			break;

		default : throw std::runtime_error("OutputHandlerCsv::add() : unhandled ColType '"+colType2string[*col]+"'");
		}
	}
	out <<'\n';
}

////////////////////////////////////////////////////////////////////////


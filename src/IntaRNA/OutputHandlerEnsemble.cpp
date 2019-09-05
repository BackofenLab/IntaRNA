
#include "IntaRNA/OutputHandlerEnsemble.h"

#include <sstream>
#include <iomanip>
#include <numeric>
#include <algorithm>

#if INTARNA_MULITHREADING
	#include <omp.h>
#endif

namespace IntaRNA {

////////////////////////////////////////////////////////////////////////////

OutputHandlerEnsemble::
OutputHandlerEnsemble(
		const OutputConstraint & outConstraint,
		std::ostream & out,
		const InteractionEnergy & energy
		)
 :	OutputHandler(outConstraint)
	, out(out)
	, energy(energy)
{
}

////////////////////////////////////////////////////////////////////////////

OutputHandlerEnsemble::
~OutputHandlerEnsemble()
{

	// ensure outputs do not intervene
	std::stringstream outTmp;
	outTmp
		<<"id1 " <<energy.getAccessibility1().getSequence().getId() <<'\n'
		<<"id2 " <<energy.getAccessibility2().getSequence().getId() <<'\n'
		<<std::setprecision(2)
		<<std::fixed // avoid scientific output format
		<<"Zall " <<getZ() <<'\n'
		<<"Eall " <<E_2_Ekcal(energy.getE(getZ())) <<'\n'
		;
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_outputStreamUpdate)
#endif
	{
		out <<outTmp.str();
	} // omp critical(intarna_omp_outputStreamUpdate)

	out.flush();
}

////////////////////////////////////////////////////////////////////////////

void
OutputHandlerEnsemble::
add( const Interaction & i )
{

	// count the interaction
	reportedInteractions++;

	// no individual interaction output
}


////////////////////////////////////////////////////////////////////////////

} // namespace

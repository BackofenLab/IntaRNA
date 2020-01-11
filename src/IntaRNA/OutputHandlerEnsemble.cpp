
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
		<<"RT "<<energy.getRT() <<'\n'
		<<std::setprecision(2)
		<<std::fixed // avoid scientific output format
//		<<"Zall " <<getZ() <<'\n'
//		<<"Zall1 " <<(E_isINF(energy.getEall1())? 0 : energy.getBoltzmannWeight(energy.getEall1())) <<'\n'
//		<<"Zall2 " <<(E_isINF(energy.getEall2())? 0 : energy.getBoltzmannWeight(energy.getEall2())) <<'\n'
		<<"Eall " <<(Z_equal(getZ(),Z_type(0)) ? 0 : E_2_Ekcal(energy.getE(getZ()))) <<'\n'
		<<"Eall1 " <<(E_isINF(energy.getEall1())? 0 : E_2_Ekcal(energy.getEall1())) <<'\n'
		<<"Eall2 " <<(E_isINF(energy.getEall2())? 0 : E_2_Ekcal(energy.getEall2())) <<'\n'
		<<"EallTotal " <<(Z_equal(getZ(),Z_type(0))||E_isINF(energy.getEall1())||E_isINF(energy.getEall2())? 0 : E_2_Ekcal(energy.getE(getZ())+energy.getEall1()+energy.getEall2())) <<'\n'
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


#include "IntaRNA/InteractionEnergy.h"

#include <cmath>

#include <iostream>

namespace IntaRNA {

////////////////////////////////////////////////////////////////////////////

InteractionEnergy::
EnergyContributions
InteractionEnergy::
getE_contributions( const Interaction & interaction ) const
{

	// temporary access to range indices
	const size_t i1 = getIndex1(*interaction.basePairs.begin());
	const size_t i2 = getIndex2(*interaction.basePairs.begin());
	const size_t j1 = getIndex1(*interaction.basePairs.rbegin());
	const size_t j2 = getIndex2(*interaction.basePairs.rbegin());

	// fill contribution data structure
	EnergyContributions contr;
	contr.init = getE_init();
	contr.ED1 = getED1( i1, j1 );
	contr.ED2 = getED2( i2, j2 );
	contr.dangleLeft =  energyWithDangles ? Z_2_E((E_2_Z(getE_danglingLeft( i1, i2 ))*getPr_danglingLeft(i1,j1,i2,j2))) : E_type(0);
	contr.dangleRight = energyWithDangles ? Z_2_E((E_2_Z(getE_danglingRight( j1, j2 ))*getPr_danglingRight(i1,j1,i2,j2))) : E_type(0);
	contr.endLeft = getE_endLeft( i1, i2 );
	contr.endRight = getE_endRight( j1, j2 );
	contr.energyAdd = getEnergyAdd();
	// compute loop energy
	contr.loops = interaction.energy
					- contr.init
					- contr.ED1
					- contr.ED2
					- contr.dangleLeft
					- contr.dangleRight
					- contr.endLeft
					- contr.endRight
					- contr.energyAdd
					;

	// final contribution distribution
	return contr;
}

////////////////////////////////////////////////////////////////////////////


} // namespace

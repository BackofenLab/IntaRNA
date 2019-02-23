#include "IntaRNA/HelixHandler.h"
#include "IntaRNA/HelixHandlerStackingOnly.h"
#include "IntaRNA/HelixHandlerUnpaired.h"

namespace IntaRNA {

HelixHandler* HelixHandler::getHelixHandler(const InteractionEnergy &energy,
											const HelixConstraint &helixConstraint,
											SeedHandler * const seedHandler) {

	if (helixConstraint.getMaxIL() == 0) {
		return new HelixHandlerStackingOnly(energy, helixConstraint, seedHandler);
	} else {
		return new HelixHandlerUnpaired(energy, helixConstraint, seedHandler);
	}
}

}
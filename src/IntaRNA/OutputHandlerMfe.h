#ifndef INTARNA_OUTPUTHANDLERMFE_H_
#define INTARNA_OUTPUTHANDLERMFE_H_

#include "IntaRNA/InteractionEnergy.h"
#include "IntaRNA/OutputHandler.h"
#include "IntaRNA/SeedHandler.h"

#include <memory>

namespace IntaRNA {

/**
 * Final-output adapter for interactions selected by ensemble energy.
 *
 * The predictor-provided energy is retained as Interaction::ensembleEnergy.
 * For each interaction that reaches this adapter, the minimum-free-energy
 * structure with the same boundaries and seed requirements is traced and its
 * energy is stored in Interaction::energy before forwarding. This keeps the
 * additional MFE computation proportional to the final reported result count.
 */
class OutputHandlerMfe : public OutputHandler {
public:
	/**
	 * Construction takes ownership of both supplied handlers.
	 * @param output final output handler receiving representative structures
	 * @param energy interaction energy model used for the ensemble prediction
	 * @param seedHandler optional seed handler; NULL denotes seed-free prediction
	 */
	OutputHandlerMfe(OutputHandler * output, const InteractionEnergy & energy,
			SeedHandler * seedHandler = NULL);

	//! Destruction of the owned output and seed handlers.
	virtual ~OutputHandlerMfe() override;

	/**
	 * Trace and forward one representative MFE interaction.
	 * @param interaction ensemble-selected interaction boundaries
	 */
	void add(const Interaction & interaction) override;

	//! @return number of interactions forwarded by the wrapped output handler
	size_t reported() const override;

	//! Forward partition totals to the wrapped output handler.
	void incrementZ(Z_type subZ) override;

	//! @return partition total of the wrapped output handler
	Z_type getZ() const override;

	//! Forward exact-partition presentation policy.
	void setExactPartition(bool exact) override;

private:
	//! final format-specific output handler
	std::unique_ptr<OutputHandler> output;
	//! unshifted energy model shared with the prediction
	const InteractionEnergy & energy;
	//! optional seed constraint used for fixed-boundary traceback
	std::unique_ptr<SeedHandler> seedHandler;

	/** Fill interaction with the representative structure and MFE. */
	void traceMfe(Interaction & interaction);
};

} // namespace IntaRNA

#endif /* INTARNA_OUTPUTHANDLERMFE_H_ */

#ifndef INTARNA_PREDICTOREVALONLY_H_
#define INTARNA_PREDICTOREVALONLY_H_

#include "IntaRNA/Predictor.h"

#include <string>
#include <vector>

namespace IntaRNA {

/**
 * Evaluates predefined, nested intermolecular base pairs without searching.
 * Prediction ranges and output filters are ignored. Identical structures are
 * counted once; Zall and trackers describe only the supplied structures.
 */
class PredictorEvalOnly : public Predictor {
public:

	/**
	 * Copies and validates the structures, discarding input energies and seeds.
	 * @param energy the energy model; it and its sequences must outlive this object
	 * @param output the reporting destination, which must outlive this object
	 * @param predTracker optional tracker owned by this predictor
	 * @param interactions nonempty list with ascending target/descending query
	 *        base-pair indices in the original, zero-based sequence coordinates
	 * @throws std::invalid_argument for empty, incompatible or invalid structures
	 */
	PredictorEvalOnly( const InteractionEnergy & energy, OutputHandler & output,
			PredictionTracker * predTracker, const std::vector<Interaction> & interactions );

	/**
	 * Evaluates all structures using initiation, loop, accessibility, dangling-end,
	 * terminal-pair and energy-shift contributions from the selected model.
	 * @param r1 ignored; every supplied interaction is evaluated
	 * @param r2 ignored; every supplied interaction is evaluated
	 * @throws std::runtime_error if the model cannot assign a finite energy
	 */
	void predict( const IndexRange & r1 = IndexRange(0,RnaSequence::lastPos),
			const IndexRange & r2 = IndexRange(0,RnaSequence::lastPos) ) override;

	/**
	 * Parses a colon-separated list of hybridDB (start1dotbar1&start2dotbar2)
	 * encodings, including full-length encodings with flanking dots. Pairing bars
	 * are matched antiparallel. Starts use each sequence's input/output indexing.
	 * @param encoding one or more nonempty encodings, each containing a base pair
	 * @param target the first sequence; must outlive the returned interactions
	 * @param query the second sequence in its original 5'-3' orientation
	 * @return validated interactions referencing the provided sequences
	 * @throws std::invalid_argument for malformed, out-of-bounds, unbalanced or
	 *         non-complementary structures
	 */
	static std::vector<Interaction> parseInteractions( const std::string & encoding,
			const RnaSequence & target, const RnaSequence & query );

protected:
	//! reset the partition function before each evaluation
	void initOptima() override;
	//! record an evaluated energy for the partition function and tracker
	void updateOptima( const size_t i1, const size_t j1,
			const size_t i2, const size_t j2, const E_type energy,
			const bool isHybridE, const bool incrementZall ) override;
	//! forward the partition function and all evaluated structures to output
	void reportOptima() override;

private:
	//! unique structures referencing the energy model's original sequences
	std::vector<Interaction> interactions;
};

} // namespace IntaRNA

#endif

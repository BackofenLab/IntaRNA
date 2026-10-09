#ifndef INTARNA_BASEPAIRPROBABILITYWRITER_H_
#define INTARNA_BASEPAIRPROBABILITYWRITER_H_
#include "IntaRNA/BasePairProbabilities.h"
#include "IntaRNA/RnaSequence.h"
#include <ostream>
#include <string>

namespace IntaRNA {
class InteractionEnergy;
/** Explicit writer for a completed actual-pair matrix, independent of trackers.
 * Whole blocks are synchronized for shared streams. No destructor emits data.
 */
class BasePairProbabilityWriter {
public:
	/** Validate and write one matrix with original 5'-to-3' nucleotide/index labels.
	 * @param out destination stream; an I/O failure throws
	 * @param result explicitly finalized raw result
	 * @param target original target sequence
	 * @param query original query sequence
	 * @param separator CSV-style separator, shared with ordinary output
	 */
	static void write(std::ostream & out,const BasePairProbabilities & result,
			const RnaSequence & target,const RnaSequence & query,const std::string & separator=";");
	/** Write to a filename or STDOUT/STDERR, using existing gzip stream support.
	 * Checks stream flush and explicit close; does not publish pending/failed data.
	 * @param filename destination (already expanded for multi-FASTA)
	 * @param result finalized raw result
	 * @param target original target sequence
	 * @param query original query sequence
	 * @param separator CSV-style separator
	 */
	static void writeFile(const std::string & filename,const BasePairProbabilities & result,
			const RnaSequence & target,const RnaSequence & query,const std::string & separator=";");
	/** Write a complete SVG dot plot with original target rows/query columns.
	 * Includes 1-nt accessibility, signed coordinate guides and a CSS color legend.
	 * Seed outlines are included if the result collected admitted seed pairs.
	 * @param out destination stream; an I/O failure throws
	 * @param result explicitly finalized raw result
	 * @param energy original, unshifted energy/accessibility model used in prediction
	 * @param targetRT positive finite energy scale used to convert target ED to Pu
	 * @param queryRT positive finite energy scale used to convert query ED to Pu;
	 * accessibility scales can differ from the interaction model's RT
	 */
	static void writeSvg(std::ostream & out,const BasePairProbabilities & result,
			const InteractionEnergy & energy,Z_type targetRT,Z_type queryRT);
	/** Write a complete SVG to a file or STDOUT/STDERR, including gzip support.
	 * Validates the whole document before opening the destination.
	 * @param filename destination (already expanded for multi-FASTA)
	 * @param result explicitly finalized raw result
	 * @param energy original, unshifted energy/accessibility model
	 * @param targetRT positive finite target accessibility energy scale
	 * @param queryRT positive finite query accessibility energy scale
	 */
	static void writeSvgFile(const std::string & filename,const BasePairProbabilities & result,
			const InteractionEnergy & energy,Z_type targetRT,Z_type queryRT);
private:
	static std::string svgBlock(const BasePairProbabilities & result,const InteractionEnergy & energy,
			Z_type targetRT,Z_type queryRT);
	static void emitFile(const std::string & name,const std::string & data);
	static std::string block(const BasePairProbabilities & result,const RnaSequence & target,
			const RnaSequence & query,const std::string & separator);
	static void emit(std::ostream & out,const std::string & data);
};
}
#endif

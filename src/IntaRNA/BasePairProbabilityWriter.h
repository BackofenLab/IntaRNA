#ifndef INTARNA_BASEPAIRPROBABILITYWRITER_H_
#define INTARNA_BASEPAIRPROBABILITYWRITER_H_
#include "IntaRNA/BasePairProbabilities.h"
#include "IntaRNA/RnaSequence.h"
#include <ostream>
#include <string>

namespace IntaRNA {
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
private:
	static std::string block(const BasePairProbabilities & result,const RnaSequence & target,
			const RnaSequence & query,const std::string & separator);
	static void emit(std::ostream & out,const std::string & data);
};
}
#endif

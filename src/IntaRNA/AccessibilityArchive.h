#ifndef INTARNA_ACCESSIBILITYARCHIVE_H_
#define INTARNA_ACCESSIBILITYARCHIVE_H_

#include "IntaRNA/Accessibility.h"
#include "IntaRNA/Matrix.h"

#include <boost/serialization/access.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/split_member.hpp>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace IntaRNA {

/**
 * Serialization view of the logical accessibility matrix, independent of its
 * physical storage. Version 1 stores exact internal ED integers, in rows of
 * increasing start position and interval length, including maxLength+1 for
 * dangling ends. Only valid cells are stored. Matrix-backed saves and full-band
 * loads use row views directly; generic saves and discarded input tails need
 * O(maxLength) scratch memory. The version 1 archive layout is unchanged.
 *
 * The input sequence and requested band bound allocations on load. Native
 * Boost binary archives require a compatible architecture and Boost version.
 * This is an internal file-format helper, not an installed public interface.
 */
class AccessibilityArchive {
public:
	typedef UpperBandedMatrix<E_type> EdMatrix;

	/**
	 * Create a read-only serialization view.
	 * @param source accessibility data, alive until serialization finishes
	 * @param matrix optional non-owning view of the source's unconstrained ED
	 * values; nullptr or non-empty constraints select the generic getED() path
	 */
	explicit AccessibilityArchive( const Accessibility & source, const EdMatrix * matrix = nullptr );

	/** Create a loading view; target is resized to the requested available band. */
	AccessibilityArchive( const RnaSequence & sequence, size_t maxLength, EdMatrix & target );

	/** Return the interaction length retained after loading. */
	size_t getMaxLength() const;

private:
	friend class boost::serialization::access;
	const RnaSequence & sequence;
	size_t maxLength;
	const Accessibility * source;
	//! Optional matrix storage for direct output; never owned or modified.
	const EdMatrix * sourceMatrix;
	EdMatrix * target;

	template<class Archive> void save( Archive & archive, unsigned int version ) const;
	template<class Archive> void load( Archive & archive, unsigned int version );
	BOOST_SERIALIZATION_SPLIT_MEMBER()
};

inline
AccessibilityArchive::AccessibilityArchive( const Accessibility & source, const EdMatrix * matrix )
	: sequence(source.getSequence()), maxLength(source.getMaxLength()), source(&source)
	, sourceMatrix(source.getAccConstraint().isEmpty() ? matrix : nullptr), target(nullptr)
{}

inline
AccessibilityArchive::AccessibilityArchive( const RnaSequence & sequence, size_t maxLength, EdMatrix & target )
	: sequence(sequence), maxLength(maxLength), source(nullptr), sourceMatrix(nullptr), target(&target)
{}

inline size_t
AccessibilityArchive::getMaxLength() const
{
	return maxLength;
}

template<class Archive> void
AccessibilityArchive::save( Archive & archive, unsigned int ) const
{
	const std::uint32_t magic = 0x49414343; // IACC: IntaRNA accessibility
	const std::uint32_t formatVersion = 1;
	const std::uint64_t length = sequence.size(), band = maxLength;
	const E_type infinity = Accessibility::ED_UPPER_BOUND;
	archive & magic & formatVersion & length & band & infinity;
	archive & boost::serialization::make_array(sequence.asString().data(), sequence.size());
	// Avoid maxLength+1 overflow when the entire sequence is covered.
	const size_t width = maxLength < sequence.size() ? maxLength+1 : sequence.size();
	if (sourceMatrix && (sourceMatrix->size1() != sequence.size() || sourceMatrix->size2() != sequence.size()))
		throw std::runtime_error("Accessibility archive: source matrix shape mismatch");
	std::vector<E_type> scratch(sourceMatrix ? 0 : width);
	for (size_t i = 0; i < sequence.size(); ++i) {
		const size_t count = std::min(width, sequence.size()-i);
		std::span<const E_type> row;
		if (sourceMatrix) {
			row = sourceMatrix->row(i);
			if (row.size() < count)
				throw std::runtime_error("Accessibility archive: source matrix band too narrow");
			row = row.first(count);
		} else {
			for (size_t k = 0; k < count; ++k) scratch[k] = source->getED(i, i+k);
			row = std::span<const E_type>(scratch.data(), count);
		}
		for (const E_type value : row)
			if (value < 0 || value > infinity)
				throw std::runtime_error("Accessibility archive: invalid ED value");
		archive & boost::serialization::make_array(row.data(), row.size());
	}
}

template<class Archive> void
AccessibilityArchive::load( Archive & archive, unsigned int )
{
	std::uint32_t magic, formatVersion;
	std::uint64_t length, band;
	E_type infinity;
	archive & magic & formatVersion & length & band & infinity;
	if (magic != 0x49414343 || formatVersion != 1)
		throw std::runtime_error("Accessibility archive: invalid signature or unsupported format version");
	if (length == 0 || length != sequence.size() || band > length)
		throw std::runtime_error("Accessibility archive: sequence length or matrix band mismatch");
	if (infinity != Accessibility::ED_UPPER_BOUND)
		throw std::runtime_error("Accessibility archive: incompatible energy representation");
	std::string storedSequence(sequence.size(), '\0');
	archive & boost::serialization::make_array(storedSequence.data(), storedSequence.size());
	if (storedSequence != sequence.asString())
		throw std::runtime_error("Accessibility archive: sequence mismatch");

	const size_t storedMaxLength = static_cast<size_t>(band);
	maxLength = std::min(maxLength, storedMaxLength);
	const size_t storedWidth = storedMaxLength < sequence.size() ? storedMaxLength+1 : sequence.size();
	const size_t width = maxLength < sequence.size() ? maxLength+1 : sequence.size();
	if (sequence.size() > std::numeric_limits<size_t>::max() / width / sizeof(E_type))
		throw std::runtime_error("Accessibility archive: matrix size overflow");
	target->resize(sequence.size(), sequence.size(), 0, width-1, false);
	// Read retained cells straight into their final storage. Only the discarded
	// suffix of a wider input band needs scratch space, and it is still validated.
	std::vector<E_type> discarded(storedWidth-width);
	for (size_t i = 0; i < sequence.size(); ++i) {
		const size_t count = std::min(storedWidth, sequence.size()-i);
		auto row = target->row(i);
		archive & boost::serialization::make_array(row.data(), row.size());
		const size_t tailSize = count-row.size();
		if (tailSize) archive & boost::serialization::make_array(discarded.data(), tailSize);
		for (const auto values : {std::span<const E_type>(row), std::span<const E_type>(discarded.data(), tailSize)})
			for (const E_type value : values)
				if (value < 0 || value > infinity)
					throw std::runtime_error("Accessibility archive: invalid ED value");
	}
}

} // namespace IntaRNA
#endif

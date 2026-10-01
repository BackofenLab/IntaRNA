#include "catch.hpp"

#include "IntaRNA/AccessibilityFromStream.h"
#include "IntaRNA/AccessibilityBasePair.h"
#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/AccessibilityVrna.h"
#include "IntaRNA/ReverseAccessibility.h"
#include <cstdint>
#include <cstring>
#include <sstream>

using namespace IntaRNA;

namespace {

// Count interface lookups without changing the underlying ED semantics.
template<class Base> class CountedAccessibility : public Base {
public:
	using Base::Base;
	mutable size_t lookups = 0;
	E_type getED(size_t from, size_t to) const override;
};

template<class Base> E_type
CountedAccessibility<Base>::getED(size_t from, size_t to) const
{
	++lookups;
	return Base::getED(from, to);
}

void checkBinaryRoundTrip( const Accessibility & source, size_t requestedLength )
{
	std::stringstream bytes;
	source.writeBinary(bytes);
	AccessibilityFromStream loaded(source.getSequence(), requestedLength, nullptr,
		bytes, AccessibilityFromStream::IntaRNA_Binary, 0.0);
	const size_t expectedLength = std::min(source.getMaxLength(),
		requestedLength == 0 ? source.getSequence().size() : requestedLength);
	REQUIRE(loaded.getMaxLength() == expectedLength);
	for (size_t i = 0; i < source.getSequence().size(); ++i) {
		for (size_t j = i; j < source.getSequence().size(); ++j) {
			const E_type expected = j-i <= expectedLength ? source.getED(i,j) : Accessibility::ED_UPPER_BOUND;
			REQUIRE(loaded.getED(i,j) == expected);
		}
	}
	// Re-exporting loaded matrices must retain the dangling-end band too.
	std::stringstream second;
	loaded.writeBinary(second);
	AccessibilityFromStream reloaded(source.getSequence(), 0, nullptr,
		second, AccessibilityFromStream::IntaRNA_Binary, 10.0);
	REQUIRE(reloaded.getMaxLength() == loaded.getMaxLength());
	for (size_t i = 0; i < source.getSequence().size(); ++i)
		for (size_t j = i; j < source.getSequence().size(); ++j)
			REQUIRE(reloaded.getED(i,j) == loaded.getED(i,j));
}

// Locate the application header without depending on the Boost archive preamble.
size_t headerOffset( const std::string & data )
{
	const std::uint32_t magic = 0x49414343;
	const size_t offset = data.find(std::string(reinterpret_cast<const char *>(&magic), sizeof(magic)));
	REQUIRE(offset != std::string::npos);
	return offset;
}

template<class T> std::string replaceField( std::string data, size_t offset, T value )
{
	std::memcpy(&data.at(offset), &value, sizeof(value));
	return data;
}

} // namespace

TEST_CASE("Binary accessibility preserves exact ED matrices", "[AccessibilityBinary]")
{
#include "testEasyLoggingSetup.icc"
	RnaSequence rna("test", "GGGGAAAACCCCUAGC");
	VrnaHandler vrna(37, "Turner04", false, false);
	for (size_t length : {size_t(1), size_t(5), rna.size()}) {
		AccessibilityVrna folded(rna, length, nullptr, vrna, rna.size());
		AccessibilityBasePair basePair(rna, length, nullptr);
		AccessibilityDisabled disabled(rna, length, nullptr);
		ReverseAccessibility reversed(folded);
		for (const Accessibility * acc : {static_cast<const Accessibility *>(&folded),
			static_cast<const Accessibility *>(&basePair), static_cast<const Accessibility *>(&disabled),
			static_cast<const Accessibility *>(&reversed)}) {
			for (size_t requested : {size_t(0), size_t(1), size_t(3), rna.size()})
				checkBinaryRoundTrip(*acc, requested);
		}
	}
	SECTION("constraints and infinity are retained in the values") {
		AccessibilityConstraint constraint(rna, "xp............px", 0, "", "", "");
		AccessibilityVrna folded(rna, 5, &constraint, vrna, rna.size());
		REQUIRE(folded.getED(1,1) == Accessibility::ED_UPPER_BOUND);
		checkBinaryRoundTrip(folded, 0);
	}
	SECTION("text input with only the dangling-end column") {
		RnaSequence shortRna("short", "AC");
		std::istringstream text("#unpaired probabilities\n #i$\tl=1\n1\t0.5\n2\t1.0\n");
		AccessibilityFromStream loaded(shortRna, 0, nullptr, text,
			AccessibilityFromStream::Pu_RNAplfold_Text, 1.0);
		REQUIRE(loaded.getMaxLength() == 0);
		checkBinaryRoundTrip(loaded, 0);
	}
	SECTION("short sequence") {
		RnaSequence shortRna("short", "A");
		AccessibilityDisabled disabled(shortRna, 0, nullptr);
		checkBinaryRoundTrip(disabled, 0);
	}
}

TEST_CASE("Binary accessibility rejects invalid archives", "[AccessibilityBinary]")
{
#include "testEasyLoggingSetup.icc"
	RnaSequence rna("test", "GGGGAAAACCCC");
	AccessibilityDisabled source(rna, 5, nullptr);
	std::stringstream stream;
	source.writeBinary(stream);
	const std::string valid = stream.str();
	const size_t header = headerOffset(valid);
	auto rejects = [&](const std::string & data) {
		std::istringstream input(data);
		REQUIRE_THROWS_AS(AccessibilityFromStream(rna, 3, nullptr, input,
			AccessibilityFromStream::IntaRNA_Binary, 1.0), std::exception);
	};
	SECTION("every truncated archive") {
		for (size_t size = 0; size < valid.size(); ++size) rejects(valid.substr(0,size));
	}
	SECTION("header and matrix values") {
		rejects("#unpaired probabilities\n");
		rejects(valid + "trailing data");
		rejects(replaceField(valid, header, std::uint32_t(0)));
		rejects(replaceField(valid, header+4, std::uint32_t(2)));
		rejects(replaceField(valid, header+8, std::uint64_t(rna.size()+1)));
		rejects(replaceField(valid, header+16, std::uint64_t(0)));
		rejects(replaceField(valid, header+16, std::uint64_t(-1)));
		rejects(replaceField(valid, header+24, E_type(1)));
		rejects(replaceField(valid, valid.size()-sizeof(E_type), E_type(-1)));
		// Invalid data in the discarded suffix of the first row must still fail.
		rejects(replaceField(valid, header+28+rna.size()+5*sizeof(E_type), E_type(-1)));
		rejects(replaceField(valid, valid.size()-sizeof(E_type), Accessibility::ED_UPPER_BOUND+1));
	}
	SECTION("different sequence of the same length") {
		RnaSequence other("other", "AAAAAAAACCCC");
		std::istringstream input(valid);
		REQUIRE_THROWS_WITH(AccessibilityFromStream(other, 5, nullptr, input,
			AccessibilityFromStream::IntaRNA_Binary, 1.0), Catch::Contains("sequence mismatch"));
	}
	SECTION("output failure") {
		std::ostringstream output;
		output.setstate(std::ios::badbit);
		REQUIRE_THROWS(source.writeBinary(output));
	}
}

TEST_CASE("Stored accessibility matrices serialize without interface lookups", "[AccessibilityBinary]")
{
#include "testEasyLoggingSetup.icc"
	RnaSequence rna("test", "GGGGAAAACCCCUAGC");
	VrnaHandler vrna(37, "Turner04", false, false);
	for (size_t length : {size_t(1), size_t(5), rna.size()}) {
		CountedAccessibility<AccessibilityVrna> folded(rna, length, nullptr, vrna, rna.size());
		std::stringstream direct, generic;
		const Accessibility & polymorphic = folded;
		polymorphic.writeBinary(direct);
		REQUIRE(folded.lookups == 0);
		folded.Accessibility::writeBinary(generic);
		REQUIRE(folded.lookups > 0);
		// Same version 1 payload, independent of matrix padding and dispatch.
		REQUIRE(direct.str() == generic.str());
		for (size_t requested : {size_t(0), size_t(1), size_t(3)}) {
			std::istringstream input(direct.str());
			CountedAccessibility<AccessibilityFromStream> loaded(rna, requested, nullptr,
				input, AccessibilityFromStream::IntaRNA_Binary, 1.0);
			std::stringstream stored, gathered;
			static_cast<const Accessibility &>(loaded).writeBinary(stored);
			REQUIRE(loaded.lookups == 0);
			loaded.Accessibility::writeBinary(gathered);
			REQUIRE(loaded.lookups > 0);
			REQUIRE(stored.str() == gathered.str());
		}
	}
}

TEST_CASE("Constrained matrix output retains getED masking", "[AccessibilityBinary]")
{
#include "testEasyLoggingSetup.icc"
	// The short-sequence matrix contains zeros even at inaccessible positions.
	RnaSequence rna("short", "AC");
	AccessibilityConstraint constraint(rna, "pb", 0, "", "", "");
	VrnaHandler vrna(37, "Turner04", false, false);
	CountedAccessibility<AccessibilityVrna> folded(rna, 1, &constraint, vrna, 0);
	std::stringstream direct, generic;
	folded.writeBinary(direct);
	REQUIRE(folded.lookups > 0);
	folded.Accessibility::writeBinary(generic);
	REQUIRE(direct.str() == generic.str());
	checkBinaryRoundTrip(folded, 0);
}

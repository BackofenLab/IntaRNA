// Compare generic/direct matrix export and raw/gzip binary I/O on identical data.
// Usage: accessibility-benchmark SEQUENCE_LENGTH OUTPUT_DIRECTORY [FASTA_FILE]
#include <IntaRNA/AccessibilityVrna.h>
#include <IntaRNA/AccessibilityFromStream.h>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

INITIALIZE_EASYLOGGINGPP

namespace {
using namespace IntaRNA;
namespace bio = boost::iostreams;
using Clock = std::chrono::steady_clock;
constexpr size_t interactionLength = 100;
constexpr size_t foldingWindow = 150;

std::string makeSequence(size_t length, const char * fasta)
{
	std::string sequence;
	sequence.reserve(length);
	if (fasta) {
		std::ifstream input(fasta);
		if (!input) throw std::runtime_error("cannot open FASTA input");
		std::string line;
		bool inSequence = false;
		while (sequence.size() < length && std::getline(input, line)) {
			if (!line.empty() && line[0] == '>') {
				if (inSequence) break;
				inSequence = true;
				continue;
			}
			for (char c : line) {
				if (c != ' ' && c != '\t' && c != '\r' && sequence.size() < length)
					sequence += c;
			}
		}
		if (sequence.size() != length) throw std::runtime_error("first FASTA sequence is too short");
	} else {
		std::uint32_t state = 245;
		for (size_t i = 0; i < length; ++i) {
			state = state * 1664525u + 1013904223u;
			sequence += "ACGU"[state >> 30];
		}
	}
	return sequence;
}

void writeCache(const Accessibility & source, const std::string & path, bool compressed, bool generic)
{
	// Match production stream buffering; only the gzip filter differs.
	constexpr std::streamsize bufferSize = 512 * 1024;
	bio::filtering_ostream output;
	if (compressed) output.push(bio::gzip_compressor(), bufferSize);
	output.push(bio::file_descriptor_sink(path, std::ios::out | std::ios::binary), bufferSize);
	if (generic) source.Accessibility::writeBinary(output);
	else source.writeBinary(output);
	output.reset(); // include compressor finalization and file close in timing
}

std::unique_ptr<AccessibilityFromStream> readCache(const RnaSequence & rna, const std::string & path, bool compressed)
{
	bio::filtering_istream input;
	if (compressed) input.push(bio::gzip_decompressor());
	input.push(bio::file_descriptor_source(path, std::ios::in | std::ios::binary));
	auto result = std::make_unique<AccessibilityFromStream>(rna, interactionLength, nullptr,
		input, AccessibilityFromStream::IntaRNA_Binary, 1.0);
	input.reset();
	return result;
}

size_t countDifferences(const Accessibility & source, const Accessibility & loaded)
{
	if (loaded.getMaxLength() != source.getMaxLength()) throw std::runtime_error("interaction length mismatch");
	size_t differences = 0;
	for (size_t i = 0; i < source.getSequence().size(); ++i)
		for (size_t j = i; j < source.getSequence().size() && j-i <= source.getMaxLength(); ++j)
			if (source.getED(i,j) != loaded.getED(i,j)) ++differences;
	return differences;
}

void benchmark(size_t length, const std::filesystem::path & directory, const char * fasta)
{
	RnaSequence rna("benchmark", makeSequence(length, fasta));
	VrnaHandler vrna(37, "Turner04", false, false);
	AccessibilityVrna folded(rna, interactionLength, nullptr, vrna, foldingWindow);
	const auto seed = (directory / "seed.bin").string();
	writeCache(folded, seed, false, false);
	auto reloaded = readCache(rna, seed, false);
	if (countDifferences(folded, *reloaded)) throw std::runtime_error("initial reload differs");
	std::filesystem::remove(seed);

	// Rotate all eight combinations so one format/method is not always first.
	std::cout << "trial,producer,serialization,format,write_seconds,read_seconds,bytes,different_cells\n";
	for (size_t trial = 0; trial < 3; ++trial) {
		for (size_t step = 0; step < 8; ++step) {
			const size_t variant = (step + 3*trial) % 8;
			const bool streamSource = variant & 4, generic = variant & 2, compressed = variant & 1;
			const Accessibility & source = streamSource ? static_cast<const Accessibility &>(*reloaded) : folded;
			const std::string producer = streamSource ? "stream" : "vrna";
			const auto path = (directory / (producer + (compressed ? ".agz" : ".bin"))).string();
			const auto start = Clock::now();
			writeCache(source, path, compressed, generic);
			const auto wrote = Clock::now();
			auto loaded = readCache(rna, path, compressed);
			const auto read = Clock::now();
			const auto differences = countDifferences(folded, *loaded);
			if (differences) throw std::runtime_error("binary ED values differ");
			std::cout << trial << ',' << producer << ',' << (generic ? "generic" : "direct") << ','
				<< (compressed ? "gzip" : "raw") << ','
				<< std::chrono::duration<double>(wrote-start).count() << ','
				<< std::chrono::duration<double>(read-wrote).count() << ','
				<< std::filesystem::file_size(path) << ',' << differences << std::endl;
		}
	}
}
} // namespace

int main(int argc, char **argv)
{
	if (argc < 3 || argc > 4) {
		std::cerr << "Usage: accessibility-benchmark SEQUENCE_LENGTH OUTPUT_DIRECTORY [FASTA_FILE]\n";
		return 2;
	}
	el::Loggers::reconfigureAllLoggers(el::ConfigurationType::Enabled, "false");
	try {
		const size_t length = std::stoull(argv[1]);
		if (length < foldingWindow) throw std::runtime_error("use at least 150 bases");
		benchmark(length, argv[2], argc == 4 ? argv[3] : nullptr);
	} catch (const std::exception & error) {
		std::cerr << error.what() << '\n';
		return 1;
	}
}

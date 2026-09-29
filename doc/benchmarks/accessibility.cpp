// Compare compressed text and binary accessibility I/O after one ViennaRNA fold.
// Usage: accessibility-benchmark SEQUENCE_LENGTH EXISTING_OUTPUT_DIRECTORY
#include <IntaRNA/AccessibilityVrna.h>
#include <IntaRNA/AccessibilityFromStream.h>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

INITIALIZE_EASYLOGGINGPP

int main(int argc, char **argv) {
    if (argc != 3) return 2;
    el::Loggers::reconfigureAllLoggers(el::ConfigurationType::Enabled, "false");
    const size_t n = std::stoull(argv[1]);
    std::string sequence;
    sequence.reserve(n);
    unsigned state = 245;
    for (size_t i=0;i<n;++i) {
        state = state * 1664525u + 1013904223u;
        sequence += "ACGU"[state >> 30];
    }
    IntaRNA::RnaSequence rna("benchmark", sequence);
    IntaRNA::VrnaHandler vrna(37, "Turner04", false, false);
    IntaRNA::AccessibilityVrna source(rna, 100, nullptr, vrna, 150);
    using Clock = std::chrono::steady_clock;
    for (int run=0;run<3;++run) {
        for (bool binary : {false, true}) {
            const auto path = std::string(argv[2]) + (binary ? "/cache.agz" : "/cache.txt.gz");
            auto start=Clock::now();
            auto *out=IntaRNA::newOutputStream(path);
            if (binary) source.writeBinary(*out); else source.writeRNAplfold_ED_text(*out);
            IntaRNA::deleteOutputStream(out);
            auto wrote=Clock::now();
            auto *in=IntaRNA::newInputStream(path);
            IntaRNA::AccessibilityFromStream read(rna,100,nullptr,*in,
                binary ? IntaRNA::AccessibilityFromStream::IntaRNA_Binary : IntaRNA::AccessibilityFromStream::ED_RNAplfold_Text, 1.0);
            IntaRNA::deleteInputStream(in);
            auto loaded=Clock::now();
            if (read.getMaxLength()!=source.getMaxLength()) return 3;
            size_t differences = 0;
            IntaRNA::E_type maxDifference = 0;
            for(size_t i=0;i<n;++i) for(size_t j=i;j<n && j<=i+100;++j) {
                const auto difference = std::abs(read.getED(i,j)-source.getED(i,j));
                if (difference != 0) ++differences;
                maxDifference = std::max(maxDifference, difference);
            }
            if (binary && differences != 0) return 4;
            std::cout << run << ',' << (binary?"agz":"text.gz") << ','
                << std::chrono::duration<double>(wrote-start).count() << ','
                << std::chrono::duration<double>(loaded-wrote).count() << ','
                << std::filesystem::file_size(path) << ',' << differences << ',' << maxDifference << std::endl;
        }
    }
}

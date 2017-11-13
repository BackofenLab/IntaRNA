
#include "catch.hpp"

#undef NDEBUG

#include <vector>
#include <utility>
#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include "IntaRNA/NussinovHandler.h"

using namespace IntaRNA;


TEST_CASE("NussinovHandler", "[NussinovHandler]") {

	// setup easylogging++ stuff if not already done
	#include "testEasyLoggingSetup.icc"

  SECTION("Dot bracket representation") {

    const std::vector<std::pair<std::string, size_t>> seqs = {{"gguccacguccaa", 3}, {"ggugcccg", 2}, {"gguccacguccaagguc", 6}};
    for (const std::pair<std::string, size_t> &seq_nuss : seqs) {
      const std::string &seq = seq_nuss.first;
      const size_t nussinovValue = seq_nuss.second;
      RnaSequence rna("test", seq);
      size_t minLoopLen = 2;
      std::string dotRep = NussinovHandler::dotBracket(0, seq.size() - 1, rna, minLoopLen);
      std::vector<size_t> open;
      size_t nuss = 0;
      for (size_t i = 0; i < dotRep.size(); ++i) {
        if (dotRep[i] == '(') {
          open.push_back(i);
          ++nuss;
        }
        if (dotRep[i] == ')') {
          REQUIRE_FALSE(open.empty());
          size_t j = open.back();
          REQUIRE(i > minLoopLen + j);
          REQUIRE(RnaSequence::areComplementary(rna, rna, j, i));
          open.pop_back();
        }
      }
      REQUIRE(nuss == nussinovValue);
    }
    //                       0123456789ABCDEFG
    //                       _______[_______]_
    const std::string seq = "gguccacguccaagguc";
    REQUIRE(seq.size() == 17);
    RnaSequence rna("test", seq);
    size_t minLoopLen = 2;
    std::string dotRep = NussinovHandler::dotBracket(7, 15, rna, minLoopLen);
    std::vector<size_t> open;
    size_t nuss = 0;
    for (size_t i = 0; i < dotRep.size(); ++i) {
      if (dotRep[i] == '(') {
        open.push_back(i);
        ++nuss;
      }
      if (dotRep[i] == ')') {
        REQUIRE_FALSE(open.empty());
        size_t j = open.back();
        REQUIRE(i > minLoopLen + j);
        REQUIRE(RnaSequence::areComplementary(rna, rna, j + 7, i + 7));
        open.pop_back();
      }
    }
    REQUIRE(nuss == 3u);

  }
}

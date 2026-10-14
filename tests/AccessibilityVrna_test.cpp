
#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/AccessibilityVrna.h"
#include "IntaRNA/InteractionEnergyVrna.h"
#include "IntaRNA/ReverseAccessibility.h"

#include <exception>
#include <memory>
#include <thread>
#include <vector>

#if INTARNA_MULITHREADING
#include <omp.h>
#endif

using namespace IntaRNA;



TEST_CASE("AccessibilityVrna", "[AccessibilityVrna]") {

	// setup easylogging++ stuff if not already done
	#include "testEasyLoggingSetup.icc"


  SECTION("ED test short") {

	const std::string seq = "GG";
	RnaSequence rna("test", seq);
	VrnaHandler vrnaHandler(37,"Turner04",false,false);
	AccessibilityVrna acc(rna, 2, NULL,vrnaHandler,2);

	REQUIRE( E_equal( acc.getED(0, 0), 0 ) );
	REQUIRE( E_equal( acc.getED(1, 1), 0 ) );
	REQUIRE( E_equal( acc.getED(0, 1), 0 ) );

  }

	SECTION("C++ string boundary preserves ViennaRNA constraints") {

		const std::string seq = "GGGGAAAACCCC";
		RnaSequence rna("test", seq);
		AccessibilityConstraint constraint(rna, "xp........px", 0, "", "", "");
		VrnaHandler vrnaHandler(37, "Turner04", false, false);
		AccessibilityVrna acc(rna, 1, &constraint, vrnaHandler, rna.size());

		// 'x' is passed at both C-string boundaries and forces unpaired bases.
		REQUIRE( E_equal(acc.getED(0, 0), 0) );
		REQUIRE( E_equal(acc.getED(rna.size()-1, rna.size()-1), 0) );
		// 'p' is passed at both C-string boundaries and forces paired bases.
		REQUIRE( E_equal(acc.getED(1, 1), Accessibility::ED_UPPER_BOUND) );
		REQUIRE( E_equal(acc.getED(rna.size()-2, rna.size()-2), Accessibility::ED_UPPER_BOUND) );
		
		// check maximal length
		REQUIRE( acc.getMaxLength() == 1 );
		// check for ED values of length+1 regions required for dangling-end probability computation
		REQUIRE( E_equal(acc.getED(0, 1), Accessibility::ED_UPPER_BOUND) ); // not accessible due to 'p' at position 1
		//std::cerr << "acc.getED(2,3) = " << acc.getED(2,3) << std::endl;
		REQUIRE( E_equal(acc.getED(2,3), 610) ); // should be accessible for dangling-end probability computation of position 2 or 3

		// Exercise the same x|........|x bytes through computeES() and
		// computeIntraEall(), including the reversed accessibility path.
		AccessibilityDisabled disabledAcc(rna, rna.size(), &constraint);
		ReverseAccessibility disabledAccReversed(disabledAcc);
		InteractionEnergyVrna energy(
				disabledAcc, disabledAccReversed, vrnaHandler, 16, 16, true);
		REQUIRE( E_isNotINF(energy.getEall1()) );
		REQUIRE( E_isNotINF(energy.getEall2()) );
	}


	SECTION("InteractionEnergyVrna respects accessibility base-pair span") {

		const std::string seq = "GGGGAAAACCCC";
		RnaSequence rna("test", seq);
		AccessibilityConstraint unrestrictedConstraint(rna.size(), 0, "", "", "");
		AccessibilityConstraint shortSpanConstraint(rna.size(), 3, "", "", "");
		AccessibilityDisabled unrestrictedAcc(rna, rna.size(), &unrestrictedConstraint);
		AccessibilityDisabled shortSpanAcc(rna, rna.size(), &shortSpanConstraint);
		ReverseAccessibility unrestrictedAccReversed(unrestrictedAcc);
		ReverseAccessibility shortSpanAccReversed(shortSpanAcc);
		VrnaHandler vrnaHandler(37, "Turner04", false, false);

		InteractionEnergyVrna unrestrictedEnergy(
				unrestrictedAcc, unrestrictedAccReversed, vrnaHandler, 16, 16, true);
		InteractionEnergyVrna shortSpanEnergy(
				shortSpanAcc, shortSpanAccReversed, vrnaHandler, 16, 16, true);

		REQUIRE( E_isNotINF(unrestrictedEnergy.getES1(0, rna.size()-1)) );
		REQUIRE( E_isINF(shortSpanEnergy.getES1(0, rna.size()-1)) );

		// Exercise both normal computeIntraEall() paths. Their ViennaRNA-owned
		// resources are released by the scoped owners in InteractionEnergyVrna.
		REQUIRE( E_isNotINF(unrestrictedEnergy.getEall1()) );
		REQUIRE( E_isNotINF(unrestrictedEnergy.getEall2()) );
	}
}

TEST_CASE("Concurrent ViennaRNA accessibility and energy setup", "[AccessibilityVrna][threads]") {

	#include "testEasyLoggingSetup.icc"

	// Configure models before starting workers, as the CLI does. Different
	// temperatures and window/span sizes expose cross-compound cache reuse.
	constexpr size_t caseCount = 12;
	constexpr size_t repeats = 4;
	constexpr size_t workerCount = 4;
	std::vector<RnaSequence> sequences;
	std::vector<std::unique_ptr<VrnaHandler>> handlers;
	std::vector<std::unique_ptr<AccessibilityConstraint>> constraints;
	sequences.reserve(caseCount);
	for (size_t c = 0; c < caseCount; ++c) {
		std::string sequence = "GGGGAAAACCCCUAGC";
		for (size_t r = 0; r < c % 3; ++r) {
			sequence += "GGGGAAAACCCCUAGC";
		}
		sequences.emplace_back("thread-case-" + std::to_string(c), sequence);
		handlers.emplace_back(std::make_unique<VrnaHandler>(15.0 + 5.0*c, "Turner04", c % 2, c % 3 == 0));
		std::string structure(sequence.size(), '.');
		if (c % 2) {
			structure.front() = structure.back() = 'x';
			structure[sequence.size()/2] = 'b';
		}
		constraints.emplace_back(std::make_unique<AccessibilityConstraint>(
				sequences.back(), structure, c % 2 ? 10 : 0, "", "", ""));
	}

	// Snapshot scientific outputs, including the extra ED band for dangles,
	// and global folding paths that can overlap target accessibility in the CLI.
	auto calculate = [&](const size_t c) {
		const RnaSequence & sequence = sequences[c];
		AccessibilityVrna acc(sequence, 8, constraints[c].get(), *handlers[c],
				c % 2 ? 12 : sequence.size(), c % 2 ? 1.2 : 1.07);
		ReverseAccessibility reverse(acc);
		InteractionEnergyVrna energy(acc, reverse, *handlers[c], 16, 16, true);
		std::vector<E_type> values;
		for (size_t i = 0; i < sequence.size(); ++i) {
			for (size_t j = i; j < sequence.size(); ++j) {
				if (j-i <= 8) {
					values.push_back(acc.getED(i,j));
				}
				values.push_back(energy.getES1(i,j));
				values.push_back(energy.getES2(i,j));
			}
		}
		values.push_back(energy.getE_init());
		values.push_back(energy.getEall1());
		values.push_back(energy.getEall2());
		return values;
	};
	std::vector<std::vector<E_type>> expected(caseCount);
	for (size_t c = 0; c < caseCount; ++c) {
		expected[c] = calculate(c);
	}

	// Catch assertions stay on the calling thread; each job owns its results
	// and exception slot, so test bookkeeping adds no data races of its own.
	std::vector<std::vector<E_type>> actual(caseCount * repeats);
	std::vector<std::exception_ptr> errors(actual.size());
	auto runJob = [&](const size_t job) {
		try {
			actual[job] = calculate(job % caseCount);
		} catch (...) {
			errors[job] = std::current_exception();
		}
	};

	SECTION("native threads also protect builds without IntaRNA OpenMP") {
		std::vector<std::jthread> workers;
		for (size_t worker = 0; worker < workerCount; ++worker) {
			workers.emplace_back([&, worker] {
				for (size_t job = worker; job < actual.size(); job += workerCount) {
					runJob(job);
				}
			});
		}
		// jthread joins before the common assertions below.
	}

#if INTARNA_MULITHREADING
	SECTION("OpenMP workers match the CLI execution model") {
		#pragma omp parallel for schedule(dynamic) num_threads(workerCount)
		for (size_t job = 0; job < actual.size(); ++job) {
			runJob(job);
		}
	}
#endif

	for (size_t job = 0; job < actual.size(); ++job) {
		CAPTURE(job);
		REQUIRE(errors[job] == nullptr);
		REQUIRE(actual[job] == expected[job % caseCount]);
	}
}

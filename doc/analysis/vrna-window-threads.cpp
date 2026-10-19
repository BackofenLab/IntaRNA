/**
 * Standalone ViennaRNA diagnostic for IntaRNA issue #265.
 *
 * Usage: vrna-window-threads MODE THREADS REPETITIONS [shared-model]
 * MODE is setup-locked, window-only, unsafe, window-locked, or setup-only.
 * The unsafe modes intentionally access ViennaRNA's shared parameter caches
 * concurrently. Run them as separate diagnostic processes, never in normal CI.
 * Passing this finite stress test does not establish thread safety.
 */
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <omp.h>

extern "C" {
#include <ViennaRNA/vrna_config.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/LPfold.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/constraints/basic.h>
#include <ViennaRNA/constraints/soft.h>
#include <ViennaRNA/params/basic.h>
#include <ViennaRNA/params/io.h>
}

namespace {

constexpr unsigned int options = VRNA_OPTION_PF | VRNA_OPTION_WINDOW;
using Compound = std::unique_ptr<vrna_fold_compound_t, decltype(&vrna_fold_compound_free)>;

struct Input {
	std::string sequence;
	std::string constraint;
	vrna_md_t model;
	int window;
	bool softConstraint;
};

struct Result {
	int length;
	int window;
	int status = 0;
	int callbacks = 0;
	bool invalid = false;
	std::vector<double> parameters;
	std::vector<double> probabilities;
	std::vector<unsigned char> seen;

	Result(int length, int window)
		: length(length), window(window),
		  probabilities((length + 1) * (window + 1)), seen(probabilities.size())
	{}
};

std::vector<Input>
makeInputs(bool sharedModel)
{
	std::vector<Input> inputs;
	uint32_t rng = 0x265;
	for (int id = 0; id < 24; ++id) {
		Input input;
		const int length = 77 + id * 13;
		for (int i = 0; i < length; ++i) {
			rng ^= rng << 13;
			rng ^= rng >> 17;
			rng ^= rng << 5;
			input.sequence += "ACGU"[rng % 4];
		}
		if (id % 4 == 0) {
			input.sequence.replace(10, 20, "GGGGGGGAAAAAACCCCCCC");
		}
		input.window = std::min(length, 32 + (id % 7) * 16);
		vrna_md_set_default(&input.model);
		input.model.window_size = input.window;
		input.model.max_bp_span = id % 3 == 0 ? input.window - 8 : input.window;
		input.model.temperature = sharedModel ? 37 : 20 + (id % 5) * 8.5;
		input.model.noLP = sharedModel ? 0 : id % 2;
		input.model.noGUclosure = sharedModel ? 0 : (id / 2) % 2;
		input.model.dangles = 2;
		input.model.gquad = 0;
		input.model.compute_bpp = 1;
		if (id % 3 == 1) {
			input.constraint.assign(input.sequence.size(), '.');
			input.constraint.front() = input.constraint.back() = 'x';
		}
		input.softConstraint = id % 3 == 2;
		inputs.push_back(input);
	}
	return inputs;
}

Compound
makeCompound(const Input &input)
{
	vrna_md_t model = input.model;
	Compound compound(vrna_fold_compound(input.sequence.c_str(), &model, options),
	                  &vrna_fold_compound_free);
	if (!compound) {
		std::cerr << "Fold compound allocation failed\n";
		std::exit(3);
	}
	if (!input.constraint.empty()) {
		vrna_constraints_add(compound.get(), input.constraint.c_str(),
		                     VRNA_CONSTRAINT_DB_DEFAULT | VRNA_CONSTRAINT_DB_ENFORCE_BP);
	}
	if (input.softConstraint) {
		// Same per-position pseudoenergy storage used by unpaired SHAPE constraints.
		vrna_sc_add_up(compound.get(), 10, 0.5, options);
		vrna_sc_add_up(compound.get(), 20, -0.25, options);
	}
	// probs_window prepares again. Do the potentially cache-writing preparation
	// here so that setup-locked includes model sanitization and parameter updates.
	if (!vrna_fold_compound_prepare(compound.get(), options)) {
		std::cerr << "Fold compound preparation failed\n";
		std::exit(3);
	}
	return compound;
}

std::vector<double>
parameterSnapshot(const vrna_fold_compound_t &compound)
{
	std::vector<double> values;
	for (const vrna_md_t *model : { &compound.params->model_details,
	                              &compound.exp_params->model_details }) {
		values.insert(values.end(), { model->temperature,
			static_cast<double>(model->window_size), static_cast<double>(model->max_bp_span),
			static_cast<double>(model->min_loop_size), static_cast<double>(model->noLP),
			static_cast<double>(model->noGUclosure), static_cast<double>(model->dangles) });
	}
	// Stable fields only: parameter IDs and structure padding are deliberately
	// excluded; they need not be identical between independently created objects.
	values.push_back(compound.exp_params->kT);
	for (int i = 0; i <= NBPAIRS; ++i) {
		for (int j = 0; j <= NBPAIRS; ++j) {
			values.push_back(compound.params->stack[i][j]);
			values.push_back(compound.exp_params->expstack[i][j]);
		}
	}
	return values;
}

void
capture(FLT_OR_DBL *probabilities, int size, int j, int max,
        unsigned int type, void *data)
{
	auto &result = *static_cast<Result *>(data);
	if (!(type & VRNA_PROBS_WINDOW_UP)) {
		return;
	}
	++result.callbacks;
	if (j < 1 || j > result.length) {
		result.invalid = true;
		return;
	}
	for (int length = 1; length <= std::min({ j, size, max, result.window }); ++length) {
		const size_t index = j * (result.window + 1) + length;
		const double probability = probabilities[length];
		result.probabilities[index] = probability;
		result.seen[index] = 1;
		if (!std::isfinite(probability) || probability < 0 || probability > 1.00000001) {
			result.invalid = true;
		}
	}
}

void
calculate(vrna_fold_compound_t *compound, const Input &input, Result &result)
{
	result.status = vrna_probs_window(compound, input.window,
	                                VRNA_PROBS_WINDOW_UP, capture, &result);
}

} // namespace

int
main(int argc, char **argv)
{
	if (argc < 4 || argc > 5) {
		std::cerr << "Usage: " << argv[0]
		          << " MODE THREADS REPETITIONS [shared-model]\n"
		          << "MODE: setup-locked | window-only | unsafe | window-locked | setup-only\n";
		return 2;
	}
	const std::string mode = argv[1];
	const int threads = std::atoi(argv[2]);
	const int repetitions = std::atoi(argv[3]);
	const bool sharedModel = argc == 5 && std::string(argv[4]) == "shared-model";
	if (threads < 1 || repetitions < 1 ||
	    (argc == 5 && !sharedModel) ||
	    (mode != "setup-locked" && mode != "window-only" && mode != "unsafe" &&
	     mode != "window-locked" && mode != "setup-only")) {
		std::cerr << "Invalid arguments\n";
		return 2;
	}
	omp_set_dynamic(0);
	omp_set_num_threads(threads);
	vrna_params_load_RNA_Turner2004();
	const auto inputs = makeInputs(sharedModel);
	std::vector<Result> reference;
	for (const auto &input : inputs) {
		reference.emplace_back(input.sequence.size(), input.window);
		auto compound = makeCompound(input);
		reference.back().parameters = parameterSnapshot(*compound);
		if (mode != "setup-only") {
			calculate(compound.get(), input, reference.back());
			if (!reference.back().status || reference.back().invalid) {
				std::cerr << "Invalid serial baseline\n";
				return 3;
			}
		}
	}

	size_t compared = 0, probabilityMismatches = 0, parameterMismatches = 0, metadataMismatches = 0;
	const auto start = std::chrono::steady_clock::now();
	int completed = 0;
	for (int repetition = 0; repetition < repetitions; ++repetition) {
		std::vector<Result> observed;
		std::vector<Compound> compounds;
		for (const auto &input : inputs) {
			observed.emplace_back(input.sequence.size(), input.window);
			compounds.emplace_back(mode == "window-only" ? makeCompound(input).release() : nullptr,
			                       &vrna_fold_compound_free);
		}
		#pragma omp parallel for schedule(dynamic, 1)
		for (int id = 0; id < static_cast<int>(inputs.size()); ++id) {
			if (mode == "setup-locked") {
				#pragma omp critical(vrna_probe_setup)
				compounds[id] = makeCompound(inputs[id]);
			} else if (mode != "window-only") {
				compounds[id] = makeCompound(inputs[id]);
			}
			observed[id].parameters = parameterSnapshot(*compounds[id]);
			if (mode == "window-locked") {
				#pragma omp critical(vrna_probe_window)
				calculate(compounds[id].get(), inputs[id], observed[id]);
			} else if (mode != "setup-only") {
				calculate(compounds[id].get(), inputs[id], observed[id]);
			}
			compounds[id].reset();
		}
		++completed;
		for (size_t id = 0; id < inputs.size(); ++id) {
			const auto &expected = reference[id];
			const auto &actual = observed[id];
			if (expected.parameters != actual.parameters) {
				if (parameterMismatches < 10) {
					const auto first = std::mismatch(expected.parameters.begin(), expected.parameters.end(),
					                                 actual.parameters.begin());
					std::cerr << std::setprecision(17) << "Parameter mismatch: repetition=" << repetition
					          << " case=" << id << " field=" << first.first - expected.parameters.begin()
					          << " serial=" << *first.first << " parallel=" << *first.second << '\n';
				}
				++parameterMismatches;
			}
			if (expected.status != actual.status || actual.invalid ||
			    expected.callbacks != actual.callbacks || expected.seen != actual.seen) {
				std::cerr << "Callback metadata mismatch: repetition=" << repetition << " case=" << id << '\n';
				++metadataMismatches;
			}
			for (size_t index = 0; index < expected.probabilities.size(); ++index) {
				if (!expected.seen[index]) {
					continue;
				}
				++compared;
				if (expected.probabilities[index] != actual.probabilities[index]) {
					if (probabilityMismatches < 10) {
						std::cerr << std::setprecision(17) << "Probability mismatch: repetition=" << repetition
						          << " case=" << id << " end=" << index / (expected.window + 1)
						          << " length=" << index % (expected.window + 1)
						          << " serial=" << expected.probabilities[index]
						          << " parallel=" << actual.probabilities[index] << '\n';
					}
					++probabilityMismatches;
				}
			}
		}
		if (parameterMismatches || probabilityMismatches || metadataMismatches) {
			break;
		}
	}
	const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
	std::cout << "ViennaRNA=" << VRNA_VERSION << " mode=" << mode << " threads=" << threads
	          << " repetitions=" << completed << '/' << repetitions << " cases=" << inputs.size()
	          << " profile=" << (sharedModel ? "shared-model" : "mixed-model")
	          << " compared=" << compared << " parameter_mismatches=" << parameterMismatches
	          << " probability_mismatches=" << probabilityMismatches
	          << " metadata_mismatches=" << metadataMismatches << " elapsed_seconds=" << elapsed << '\n';
	return parameterMismatches || probabilityMismatches || metadataMismatches ? 1 : 0;
}

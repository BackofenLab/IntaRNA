#include "IntaRNA/AccessibilityVrna.h"
#include "IntaRNA/InteractionEnergyVrna.h"
#include "IntaRNA/ReverseAccessibility.h"
#include <ViennaRNA/LPfold.h>
#include <ViennaRNA/vrna_config.h>
#include <atomic>
#include <barrier>
#include <iostream>
#include <thread>
#include <vector>
INITIALIZE_EASYLOGGINGPP
static_assert(INTARNA_MULITHREADING == 0);
#ifdef _OPENMP
#error This probe must compile without OpenMP enabled in IntaRNA.
#endif
std::atomic<int> active{0}, maximum{0};
extern "C" int __real_vrna_probs_window(vrna_fold_compound_t *, int, unsigned int, vrna_probs_window_callback *, void *);
extern "C" int __wrap_vrna_probs_window(vrna_fold_compound_t *fc, int u, unsigned int o, vrna_probs_window_callback *cb, void *data) {
    const int now = ++active;
    int seen = maximum.load();
    while (seen < now && !maximum.compare_exchange_weak(seen, now)) {}
    const int result = __real_vrna_probs_window(fc, u, o, cb, data);
    --active;
    return result;
}
int main() {
    constexpr int count = 16;
    IntaRNA::VrnaHandler handler;
    std::vector<IntaRNA::RnaSequence> sequences;
    for (int i = 0; i < count; ++i) {
        std::string seq;
        for (int j = 0; j < 128 + i; ++j) seq += "GGGAAACCCUGUA"[(j + i) % 13];
        sequences.emplace_back(std::to_string(i), seq);
    }
    auto calculate = [&](int i) {
        auto &seq = sequences[i];
        IntaRNA::AccessibilityConstraint constraint(seq.size(), 55 + i % 4, "", "", "");
        IntaRNA::AccessibilityVrna acc(seq, 24, &constraint, handler, 72 + i % 4);
        std::vector<IntaRNA::E_type> result;
        for (size_t j = 0; j < seq.size(); ++j)
            for (size_t k = j; k < seq.size() && k < j + 25; ++k) result.push_back(acc.getED(j,k));
        IntaRNA::ReverseAccessibility reverse(acc);
        IntaRNA::InteractionEnergyVrna energy(acc, reverse, handler, 16, 16, true);
        result.push_back(energy.getEall1());
        result.push_back(energy.getEall2());
        result.push_back(energy.getES1(0, seq.size()-1));
        return result;
    };
    std::vector<std::vector<IntaRNA::E_type>> expected(count), actual(count);
    for (int i = 0; i < count; ++i) expected[i] = calculate(i);
    maximum = 0;
    for (int repetition = 0; repetition < 3; ++repetition) {
        std::atomic<int> next{0};
        std::barrier start(4);
        std::vector<std::thread> workers;
        for (int t = 0; t < 4; ++t) workers.emplace_back([&] {
            start.arrive_and_wait();
            for (int i; (i = next++) < count;) actual[i] = calculate(i);
        });
        for (auto &worker: workers) worker.join();
        for (int i = 0; i < count; ++i) if (actual[i] != expected[i]) return 2;
    }
    std::cout << "INTARNA_MULITHREADING=" << INTARNA_MULITHREADING
              << " threads=4 cases=16 repetitions=3 all_energy_values_equal=1 max_window_calls=" << maximum << '\n';
#ifdef VRNA_WITH_OPENMP
    std::cout << "RNAlib header OpenMP enabled (control)\n";
    return maximum > 1 ? 0 : 3;
#else
    std::cout << "RNAlib header OpenMP disabled by test override; underlying RNAlib remains OpenMP-enabled\n";
    return maximum == 1 ? 0 : 4;
#endif
}

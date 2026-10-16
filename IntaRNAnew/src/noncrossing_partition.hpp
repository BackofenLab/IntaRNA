#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace intarnanew::detail {

struct NoncrossingIntervalPartitions {
    double logPartition{};
    // Triangular interval data in begin * (n + 1) + length layout. Entries
    // that do not encode a nonempty in-range interval remain NaN.
    std::vector<double> probabilities;
};

// Computes the partition function of a noncrossing matching and the joint
// unpaired probability of every nonempty interval. Q uses half-open intervals:
//
//   Q[a,b] = u(b-1) Q[a,b-1]
//          + sum_i Q[a,i] W(i,b-1) Q[i+1,b-1].
//
// The reverse pass obtains the outside context O_R(i,j) of every pair token
// R(i,j)=W(i,j)Q[i+1,j]. For U=[a,b), structures in which U is unpaired either
// contain no pair spanning U, or have one unique innermost spanning pair:
//
//   Z_U = Q[0,a]Q[b,n]
//       + sum_{i<a,j>=b} O_R(i,j)W(i,j)Q[i+1,a]Q[b,j].
//
// Two triangular log-semiring contractions evaluate the second term for every
// U in O(n^3) total time. All arithmetic is in log space. pairLogWeight must
// return -infinity for a forbidden pair and unpairedAllowed must encode hard
// paired-position constraints.
template<class UnpairedAllowed, class PairLogWeight>
[[nodiscard]] auto computeNoncrossingIntervalPartitions(
    const std::size_t length,
    UnpairedAllowed unpairedAllowed,
    PairLogWeight pairLogWeight) -> NoncrossingIntervalPartitions {
    constexpr auto logZero = -std::numeric_limits<double>::infinity();
    const auto logAdd = [](const double left, const double right) noexcept {
        if (left == logZero) return right;
        if (right == logZero) return left;
        const auto high = std::max(left, right);
        return high + std::log1p(std::exp(std::min(left, right) - high));
    };
    const auto dimension = length + 1U;
    const auto intervalOffset = [dimension](const std::size_t begin, const std::size_t end) {
        return begin * dimension + end;
    };
    const auto pairOffset = [length](const std::size_t left, const std::size_t right) {
        return left * length + right;
    };

    std::vector<double> inside(dimension * dimension, logZero);
    std::vector<double> pairWeights(length * length, logZero);
    for (std::size_t index{}; index <= length; ++index) {
        inside[intervalOffset(index, index)] = 0.0;
    }
    for (std::size_t left{}; left < length; ++left) {
        for (std::size_t right = left + 1U; right < length; ++right) {
            pairWeights[pairOffset(left, right)] = pairLogWeight(left, right);
        }
    }

    for (std::size_t width = 1U; width <= length; ++width) {
        for (std::size_t begin{}; begin + width <= length; ++begin) {
            const auto end = begin + width;
            const auto last = end - 1U;
            auto value = unpairedAllowed(last)
                ? inside[intervalOffset(begin, last)] : logZero;
            for (std::size_t partner = begin; partner < last; ++partner) {
                const auto weight = pairWeights[pairOffset(partner, last)];
                if (weight == logZero) continue;
                const auto left = inside[intervalOffset(begin, partner)];
                const auto inner = inside[intervalOffset(partner + 1U, last)];
                if (left == logZero || inner == logZero) continue;
                value = logAdd(value, left + weight + inner);
            }
            inside[intervalOffset(begin, end)] = value;
        }
    }

    const auto logPartition = inside[intervalOffset(0U, length)];
    std::vector<double> outside(dimension * dimension, logZero);
    std::vector<double> pairOutside(length * length, logZero);
    outside[intervalOffset(0U, length)] = 0.0;
    for (std::size_t width = length; width > 0U; --width) {
        for (std::size_t begin{}; begin + width <= length; ++begin) {
            const auto end = begin + width;
            const auto last = end - 1U;
            const auto context = outside[intervalOffset(begin, end)];
            if (context == logZero) continue;
            if (unpairedAllowed(last)) {
                auto& child = outside[intervalOffset(begin, last)];
                child = logAdd(child, context);
            }
            for (std::size_t partner = begin; partner < last; ++partner) {
                const auto weight = pairWeights[pairOffset(partner, last)];
                if (weight == logZero) continue;
                const auto left = inside[intervalOffset(begin, partner)];
                const auto inner = inside[intervalOffset(partner + 1U, last)];
                if (left == logZero || inner == logZero) continue;

                auto& tokenContext = pairOutside[pairOffset(partner, last)];
                tokenContext = logAdd(tokenContext, context + left);
                auto& leftContext = outside[intervalOffset(begin, partner)];
                leftContext = logAdd(leftContext, context + weight + inner);
                auto& innerContext = outside[intervalOffset(partner + 1U, last)];
                innerContext = logAdd(innerContext, context + left + weight);
            }
        }
    }

    // rightContext[i,b] contracts every possible right endpoint j of a pair
    // spanning [a,b); the remaining left contraction depends on a.
    std::vector<double> rightContext(length * dimension, logZero);
    for (std::size_t left{}; left < length; ++left) {
        for (std::size_t boundary{}; boundary <= length; ++boundary) {
            auto value = logZero;
            for (std::size_t right = boundary; right < length; ++right) {
                const auto context = pairOutside[pairOffset(left, right)];
                const auto weight = pairWeights[pairOffset(left, right)];
                const auto suffix = inside[intervalOffset(boundary, right)];
                if (context == logZero || weight == logZero || suffix == logZero) continue;
                value = logAdd(value, context + weight + suffix);
            }
            rightContext[left * dimension + boundary] = value;
        }
    }

    NoncrossingIntervalPartitions result;
    result.logPartition = logPartition;
    result.probabilities.assign(
        length * dimension, std::numeric_limits<double>::quiet_NaN());
    for (std::size_t begin{}; begin < length; ++begin) {
        for (std::size_t end = begin + 1U; end <= length; ++end) {
            auto constrained = inside[intervalOffset(0U, begin)] +
                               inside[intervalOffset(end, length)];
            for (std::size_t left{}; left < begin; ++left) {
                const auto prefix = inside[intervalOffset(left + 1U, begin)];
                const auto context = rightContext[left * dimension + end];
                if (prefix == logZero || context == logZero) continue;
                constrained = logAdd(constrained, prefix + context);
            }
            const auto probability = constrained == logZero || logPartition == logZero
                ? 0.0
                : std::clamp(std::exp(std::min(0.0, constrained - logPartition)), 0.0, 1.0);
            result.probabilities[begin * dimension + (end - begin)] = probability;
        }
    }
    return result;
}

} // namespace intarnanew::detail

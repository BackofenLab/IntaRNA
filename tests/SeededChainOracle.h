#ifndef INTARNA_TEST_SEEDEDCHAINORACLE_H_
#define INTARNA_TEST_SEEDEDCHAINORACLE_H_

#include <algorithm>
#include <array>
#include <functional>
#include <map>
#include <vector>

// Independent exhaustive path enumeration. These test coordinates are zero-based
// with both strands increasing. Weight callbacks return long double, so expected
// values never pass through native integer energy conversion or the tested DP.
namespace seeded_oracle {
using Pair = std::array<size_t, 2>;
using Boundary = std::array<size_t, 4>;
using Chain = std::vector<Pair>;
struct Model {
	size_t n, m, span1, span2;
	bool noLP = false;
	std::vector<Chain> seeds;
	std::function<bool(Pair)> valid;
	std::function<long double(Pair)> init;
	std::function<long double(Pair, Pair)> edge, boundary;
};
struct Result {
	long double z = 0, unseeded = 0, pairCount = 0;
	std::map<Boundary, long double> boundaries;
	std::map<Pair, long double> mass;
};
inline bool stack(Pair a, Pair b) {
	return a[0]+1 == b[0] && a[1]+1 == b[1];
}
inline Result enumerate(const Model & model) {
	Result out;
	Chain chain;
	std::function<void(long double)> visit = [&](long double hybrid) {
		const Pair p = chain.front(), q = chain.back();
		bool noLP = true;
		for (size_t k=0; k<chain.size(); ++k)
			noLP &= (k>0 && stack(chain[k-1], chain[k]))
					|| (k+1<chain.size() && stack(chain[k], chain[k+1]));
		bool seeded = false;
		for (const Chain & seed : model.seeds)
			seeded |= std::search(chain.begin(), chain.end(), seed.begin(), seed.end()) != chain.end();
		if (!model.noLP || noLP) {
			const auto w = hybrid * model.boundary(p, q);
			if (seeded) {
				out.z += w;
				out.boundaries[{p[0], p[1], q[0], q[1]}] += w;
				out.pairCount += w * chain.size();
				for (Pair pair : chain) out.mass[pair] += w;
			} else out.unseeded += w;
		}
		for (size_t i=q[0]+1; i<model.n && i-p[0]<model.span1; ++i)
		for (size_t j=q[1]+1; j<model.m && j-p[1]<model.span2; ++j) {
			const Pair r{i,j};
			if (!model.valid(r)) continue;
			const auto t = model.edge(q,r);
			if (t == 0) continue;
			chain.push_back(r); visit(hybrid*t); chain.pop_back();
		}
	};
	for (size_t i=0; i<model.n; ++i) for (size_t j=0; j<model.m; ++j) {
		if (!model.valid({i,j})) continue;
		chain = {{i,j}};
		visit(model.init({i,j}));
	}
	return out;
}
}
#endif

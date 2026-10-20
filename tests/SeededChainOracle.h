#ifndef INTARNA_TEST_SEEDEDCHAINORACLE_H_
#define INTARNA_TEST_SEEDEDCHAINORACLE_H_

#include "IntaRNA/Interaction.h"
#include <algorithm>
#include <array>
#include <functional>
#include <map>
#include <vector>

// Independent exhaustive path enumeration. These test coordinates are zero-based
// with both strands increasing. Weight callbacks return long double, so expected
// values never pass through native integer energy conversion or the tested DP.
namespace seeded_oracle {
using Pair = IntaRNA::Interaction::BasePair;
using Boundary = std::array<size_t, 4>;
using Chain = std::vector<Pair>;
struct Model {
	size_t n, m, span1, span2;
	bool noLP = false;
	std::vector<Chain> seeds;
	std::function<bool(Pair)> valid;
	std::function<long double(Pair)> init;
	std::function<long double(Pair, Pair)> edge, boundary;
	// Optional frozen heuristic endpoint for the first completed declared seed.
	std::function<Pair(Pair, Pair)> selectRight;
};
struct Result {
	long double z = 0, unseeded = 0, pairCount = 0;
	std::map<Boundary, long double> boundaries;
	std::map<Pair, long double> mass, spots;
};
inline bool stack(Pair a, Pair b) {
	return a.first+1 == b.first && a.second+1 == b.second;
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
		const Chain * owner = nullptr;
		size_t firstHit = chain.size()+1;
		for (const Chain & seed : model.seeds) {
			const auto found=std::search(chain.begin(), chain.end(), seed.begin(), seed.end());
			if (found==chain.end()) continue;
			const size_t hit=size_t(found-chain.begin())+seed.size();
			if (hit<firstHit || (hit==firstHit && seed.size()<owner->size())) {
				owner=&seed;firstHit=hit;
			}
		}
		if (!model.noLP || noLP) {
			const auto w = hybrid * model.boundary(p, q);
			const bool retained=owner && (!model.selectRight || p==owner->front()
					|| q==owner->back() || q==model.selectRight(owner->front(),owner->back()));
			if (retained) {
				out.z += w;
				out.boundaries[{p.first, p.second, q.first, q.second}] += w;
				out.pairCount += w * chain.size();
				for (Pair pair : chain) out.mass[pair] += w;
				for (size_t i=p.first; i<=q.first; ++i)
				for (size_t j=p.second; j<=q.second; ++j) out.spots[{i,j}] += w;
			} else if (!owner) out.unseeded += w;
		}
		for (size_t i=q.first+1; i<model.n && i-p.first<model.span1; ++i)
		for (size_t j=q.second+1; j<model.m && j-p.second<model.span2; ++j) {
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

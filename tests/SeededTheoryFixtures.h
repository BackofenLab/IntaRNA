#ifndef INTARNA_TEST_SEEDEDTHEORYFIXTURES_H_
#define INTARNA_TEST_SEEDEDTHEORYFIXTURES_H_
#include "SeededChainOracle.h"
#include <cmath>
#include <fstream>
#include <string>
#include <stdexcept>

namespace seeded_oracle {
struct TheoryFixture {
	std::string name;
	Model model;
	Result expected;
	std::map<Boundary, long double> terminal, seeded, free;
	std::map<Pair, bool> valid;
};
inline std::vector<TheoryFixture> theoryFixtures() {
	std::ifstream file(std::string(INTARNA_TEST_DATA_DIR)+"/seeded-theory.dat");
	if (!file) throw std::runtime_error("cannot read seeded theory fixtures");
	std::vector<TheoryFixture> fixtures;
	std::string tag;
	while (file >> tag) {
		if (tag == "CASE") {
			fixtures.emplace_back();
			auto & f=fixtures.back();
			file >> f.name >> f.model.n >> f.model.m >> f.model.span1;
			f.model.span2=f.model.span1;
		} else {
			auto & f=fixtures.back();
			if (tag == "Z") file >> f.expected.z;
			else if (tag == "S") {
				Pair p; size_t len; file >> p[0] >> p[1] >> len;
				Chain s; for (size_t k=0;k<len;++k) s.push_back({p[0]+k,p[1]+k});
				f.model.seeds.push_back(s);
			} else if (tag == "B") {
				Boundary b; for (auto & i:b) file >> i;
				file >> f.terminal[b] >> f.seeded[b] >> f.free[b];
				f.expected.boundaries[b]=f.terminal[b]*f.seeded[b];
			} else if (tag == "P") {
				Pair p; file >> p[0] >> p[1] >> f.valid[p] >> f.expected.mass[p];
			}
		}
	}
	for (auto & f:fixtures) {
		f.model.valid=[valid=f.valid](Pair p){return valid.at(p);};
		f.model.init=[](Pair){return std::exp(1.L);};
		f.model.edge=[](Pair,Pair){return std::exp(1.L);};
		f.model.boundary=[b=f.terminal](Pair p,Pair q){
			auto it=b.find({p[0],p[1],q[0],q[1]});return it==b.end()?0.L:it->second;
		};
	}
	return fixtures;
}
}
#endif

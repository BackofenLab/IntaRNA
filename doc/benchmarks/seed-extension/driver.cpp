#include "IntaRNA/SeededPartitionFunction.h"
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

int main(int argc,char **argv) {
	using namespace IntaRNA;
	using K=SeededPartitionFunction;
	if(argc!=6) return 2;
	const std::string topology=argv[1],density=argv[2],dump=argv[5];
	const size_t n=std::stoul(argv[3]),span=std::min(n,size_t(24)),length=4;
	const bool pairs=std::stoi(argv[4]);
	auto valid=[&](K::Pair p) {return topology=="dense" || (p.first>p.second?p.first-p.second:p.second-p.first)<=3;};
	std::vector<StackSeedDomain::Occurrence> occurrences;
	if(density=="many") {
		for(size_t i=0;i+length<=n;++i) for(size_t j=0;j+length<=n;++j)
			if(valid({i,j})) occurrences.push_back({i,j,length});
	} else {
		for(size_t i:{n/4,n/2,3*n/4}) if(i+length<=n) occurrences.push_back({i,i,length});
	}
	const StackSeedDomain seeds(occurrences,n,n,span,span);
	struct Cell {K::Pair p,r;Z_type h,b;};
	std::vector<Cell> cells;
	size_t boundaries=0;
	Z_type checksum=0;
	K::Weights weights{
		valid,
		[](K::Pair p) {return Z_type(0.71)+Z_type((p.first+2*p.second)%7)/100;},
		[](K::Pair p,K::Pair q) {return Z_type(0.57)+Z_type((p.first+q.first+3*p.second+q.second)%11)/100;},
		[](K::Pair p,K::Pair q) {return Z_type(1)/(1+Z_type(q.first-p.first+q.second-p.second));},
		[&](K::Pair p,K::Pair r,Z_type h,Z_type b) {
			++boundaries; checksum+=h*b;
			if(dump!="-") cells.push_back({p,r,h,b});
		}
	};
	const K::Domain domain{n,n,span,span,3,3,false};
	const auto begin=std::chrono::steady_clock::now();
	const auto result=K::compute(domain,seeds,weights,pairs);
	const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
	Z_type moment=0;
	if(pairs) for(size_t i=0;i<n;++i) for(size_t j=0;j<n;++j) moment+=result.mass(i,j);
	std::cout<<std::setprecision(17)<<seconds<<'\t'<<result.z<<'\t'<<moment<<'\t'<<boundaries<<'\t'<<checksum<<'\t'<<occurrences.size()<<'\n';
	if(dump!="-") {
		std::ofstream output(dump);
		output<<std::setprecision(17)<<"z "<<result.z<<'\n';
		for(const auto &c:cells) output<<"b "<<c.p.first<<' '<<c.p.second<<' '<<c.r.first<<' '<<c.r.second<<' '<<c.h<<' '<<c.b<<'\n';
		if(pairs) for(size_t i=0;i<n;++i) for(size_t j=0;j<n;++j) output<<"p "<<i<<' '<<j<<' '<<result.mass(i,j)<<'\n';
		if(!output) return 3;
	}
}

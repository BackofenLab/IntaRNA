#include "IntaRNA/SeededPartitionFunction.h"
#include <chrono>
#include <cstdlib>
#include <new>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// Count ordinary C++ allocations only inside compute(), including its callbacks.
// All benchmark/kernel types use fundamental alignment. Seed construction and
// output serialization happen outside this scope; RSS is measured separately.
namespace {
bool countAllocations=false;
size_t allocationCalls=0,allocatedBytes=0;
}
void * operator new(std::size_t size) {
	void * memory=std::malloc(size?size:1);
	if(!memory) throw std::bad_alloc();
	if(countAllocations) { ++allocationCalls; allocatedBytes+=size; }
	return memory;
}
void * operator new[](std::size_t size) { return ::operator new(size); }
void operator delete(void * memory) noexcept { std::free(memory); }
void operator delete[](void * memory) noexcept { std::free(memory); }
void operator delete(void * memory,std::size_t) noexcept { std::free(memory); }
void operator delete[](void * memory,std::size_t) noexcept { std::free(memory); }

int main(int argc,char **argv) {
	using namespace IntaRNA;
	using K=SeededPartitionFunction;
	if(argc!=6 && argc!=7) return 2;
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
	if(argc==7 && std::string(argv[6])=="heuristic") weights.selectRight=[&](K::Pair,K::Pair q) {
		return K::Pair{std::min(n-1,q.first+2+(q.second%3)),std::min(n-1,q.second+2+(q.first%3))};
	};
	const K::Domain domain{n,n,span,span,3,3,false};
	const auto begin=std::chrono::steady_clock::now();
	countAllocations=true;
	const auto result=K::compute(domain,seeds,weights,pairs);
	countAllocations=false;
	const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
	Z_type moment=0;
	if(pairs) for(size_t i=0;i<n;++i) for(size_t j=0;j<n;++j) moment+=result.mass(i,j);
	std::cout<<std::setprecision(17)<<seconds<<'\t'<<result.z<<'\t'<<moment<<'\t'<<boundaries<<'\t'<<checksum<<'\t'<<occurrences.size()<<'\t'<<allocationCalls<<'\t'<<allocatedBytes<<'\n';
	if(dump!="-") {
		std::ofstream output(dump);
		output<<std::setprecision(17)<<"z "<<result.z<<'\n';
		for(const auto &c:cells) output<<"b "<<c.p.first<<' '<<c.p.second<<' '<<c.r.first<<' '<<c.r.second<<' '<<c.h<<' '<<c.b<<'\n';
		if(pairs) for(size_t i=0;i<n;++i) for(size_t j=0;j<n;++j) output<<"p "<<i<<' '<<j<<' '<<result.mass(i,j)<<'\n';
		if(!output) return 3;
	}
}

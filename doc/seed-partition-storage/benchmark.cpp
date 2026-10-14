#include "IntaRNA/SeededPartitionFunction.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <vector>

INITIALIZE_EASYLOGGINGPP

namespace {
using Kernel = IntaRNA::SeededPartitionFunction;
using IntaRNA::Z_type;

// SplitMix64 gives a stable mask without implementation-specific distributions.
uint64_t mix(uint64_t x) {
	x += UINT64_C(0x9e3779b97f4a7c15);
	x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
	x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
	return x ^ (x >> 31);
}

struct Options {
	size_t n = 64, length = 7;
	unsigned density = 100;
	bool multiple = false, pairs = false, noLP = false;
	std::string dump;
};

Options parse(int argc, char **argv) {
	Options o;
	for (int a=1; a<argc; ++a) {
		const std::string key(argv[a]);
		if (a+1 == argc) throw std::invalid_argument("missing value for " + key);
		const std::string value(argv[++a]);
		if (key=="--n") o.n=std::stoull(value);
		else if (key=="--length") o.length=std::stoull(value);
		else if (key=="--density") o.density=std::stoul(value);
		else if (key=="--starts" && (value=="single" || value=="multiple")) o.multiple=value=="multiple";
		else if (key=="--pairs" && (value=="0" || value=="1")) o.pairs=value=="1";
		else if (key=="--noLP" && (value=="0" || value=="1")) o.noLP=value=="1";
		else if (key=="--dump") o.dump=value;
		else throw std::invalid_argument("unknown option or invalid value: " + key + " " + value);
	}
	const size_t offset=o.multiple?8:0;
	if (!o.n || o.n>4096 || o.n<=offset || !o.length || o.length>o.n-offset || o.density>10000)
		throw std::invalid_argument("require 1<=n<=4096, 1<=length<=n-offset, and 0<=density<=10000 (basis points)");
	return o;
}
}

int main(int argc, char **argv) {
	try {
		const auto o=parse(argc,argv);
		const size_t offset=o.multiple?8:0;
		// Density is in basis points: 100, 1000, 3750, 10000. Force the
		// seed and its diagonal prefix valid; a second pair allows singleton
		// seeds under noLP (provided the input dimension is at least two).
		std::vector<unsigned char> valid(o.n*o.n);
		for (size_t q=0; q<valid.size(); ++q) valid[q]=mix(q+262)%10000<o.density;
		for (size_t k=0; k<std::min(o.n,offset+std::max(size_t(2),o.length)); ++k) valid[k*o.n+k]=1;
		const size_t validCount=std::count(valid.begin(),valid.end(),1);
		const std::vector<IntaRNA::StackSeedDomain::Occurrence> occurrences{{offset,offset,o.length}};
		const IntaRNA::StackSeedDomain seeds(occurrences,o.n,o.n,o.n,o.n);
		const Kernel::Domain domain{o.n,o.n,o.n,o.n,2,2,o.noLP};
		std::ofstream dump;
		if (!o.dump.empty()) {
			dump.open(o.dump);
			if (!dump) throw std::runtime_error("cannot open dump: " + o.dump);
			dump << std::hexfloat;
		}
		size_t callbacks=0, eligibleStarts=0, admittedStarts=0, maxCells=0, maxValid=0;
		// Report actual per-start dimensions and occupancy outside the core timer.
		seeds.forEachStart([&](size_t i,size_t j) {
			++eligibleStarts;
			if (!valid[i*o.n+j]) return;
			++admittedStarts;
			size_t count=0;
			for (size_t a=i; a<o.n; ++a) for (size_t b=j; b<o.n; ++b) count+=valid[a*o.n+b];
			maxCells=std::max(maxCells,(o.n-i)*(o.n-j));
			maxValid=std::max(maxValid,count);
		});
		// Binary fractions and bounded positive transitions avoid decimal parsing
		// differences and explosive growth while exercising stacks and loop edges.
		Kernel::Weights weights{
			[&](auto p){return bool(valid[p.first*o.n+p.second]);},
			[](auto p){return Z_type(0.5)+Z_type((p.first+3*p.second)%4)/16;},
			[](auto p,auto q){
				return q.first==p.first+1 && q.second==p.second+1
					? Z_type(0.625)+Z_type((p.first+q.second)%2)/8
					: Z_type(1+((p.first+q.second)%2))/64;
			},
			[](auto p,auto q){return Z_type(0.5)+Z_type((p.first+q.second)%4)/16;},
			[&](auto p,auto q,auto h,auto b){
				++callbacks;
				if (dump.is_open()) dump << "complete " << p.first << ' ' << p.second << ' '
					<< q.first << ' ' << q.second << ' ' << h << ' ' << b << '\n';
			}
		};
		const auto start=std::chrono::steady_clock::now();
		const auto result=Kernel::compute(domain,seeds,weights,o.pairs);
		const auto stop=std::chrono::steady_clock::now();
		struct rusage usage{};
		if (getrusage(RUSAGE_SELF,&usage)) throw std::runtime_error("getrusage failed");
		// Correctness dumps preserve every callback in order and every matrix cell,
		// including zero values. These runs are never included in timing samples.
		if (dump.is_open()) {
			dump << "z " << result.z << '\n';
			dump << "mass-size " << result.mass.storageSize() << '\n';
			if (o.pairs) for (size_t i=0; i<o.n; ++i) for (size_t j=0; j<o.n; ++j)
				dump << "mass " << i << ' ' << j << ' ' << result.mass(i,j) << '\n';
			dump.close();
			if (!dump) throw std::runtime_error("writing dump failed: " + o.dump);
		}
		const auto elapsed=std::chrono::duration<double,std::milli>(stop-start).count();
		const auto rssKiB=usage.ru_maxrss
#ifdef __APPLE__
			/1024
#endif
			;
		const size_t states=std::max(size_t(2),o.length)+2;
		// Dense h/dh plus total/dt, valid, ending and output mass; excludes
		// Matrix objects, seed maps, the input mask, allocator and process overhead.
		const size_t denseBytes=maxCells*((states+2)*sizeof(Z_type)*(o.pairs?2:1)+1+sizeof(size_t))
			+(o.pairs?o.n*o.n*sizeof(Z_type):0);
		std::cout << std::setprecision(17)
			<< "{\"n\":" << o.n << ",\"density_bp\":" << o.density
			<< ",\"valid_pairs\":" << validCount << ",\"length\":" << o.length
			<< ",\"starts\":\"" << (o.multiple?"multiple":"single")
			<< "\",\"pairs\":" << o.pairs << ",\"noLP\":" << o.noLP
			<< ",\"eligible_starts\":" << eligibleStarts << ",\"admitted_starts\":" << admittedStarts
			<< ",\"max_cells\":" << maxCells << ",\"max_valid\":" << maxValid
			<< ",\"dense_payload_bytes\":" << denseBytes << ",\"callbacks\":" << callbacks
			<< ",\"z\":\"" << std::hexfloat << result.z << std::defaultfloat
			<< "\",\"core_ms\":" << elapsed << ",\"peak_rss_kib\":" << rssKiB << "}\n";
		return 0;
	} catch (const std::exception &e) {
		std::cerr << "benchmark: " << e.what() << '\n';
		return 1;
	}
}

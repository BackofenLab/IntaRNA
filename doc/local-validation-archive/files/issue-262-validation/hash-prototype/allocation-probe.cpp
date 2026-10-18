// Count allocation requests for the hash prototype, excluding allocator overhead.
// This is a layout probe, not a timing benchmark or a measurement of RSS.
#include "IntaRNA/Interaction.h"
#include <array>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <type_traits>
#include <unordered_map>
#include <vector>

using Pair = IntaRNA::Interaction::BasePair;
using Z = IntaRNA::Z_type;
struct Counter { size_t calls = 0, bytes = 0; };
std::array<Counter, 3> requested;

template<class T>
struct CountingAllocator {
	using value_type = T;
	CountingAllocator() = default;
	template<class U> CountingAllocator(const CountingAllocator<U> &) {}
	T *allocate(size_t n) {
		T *result = std::allocator<T>{}.allocate(n);
		// libstdc++ rebinds the map allocator to nodes and bucket pointers.
		auto &counter = requested[std::is_same_v<T, Z> ? 0 : std::is_pointer_v<T> ? 1 : 2];
		++counter.calls;
		counter.bytes += n * sizeof(T);
		return result;
	}
	void deallocate(T *p, size_t n) { std::allocator<T>{}.deallocate(p, n); }
	template<class U> bool operator==(const CountingAllocator<U> &) const { return true; }
};

int main(int argc, char **argv) {
	const size_t admitted = argc > 1 ? std::stoull(argv[1]) : 10000;
	if (!admitted) return 1;
	using Values = std::vector<Z, CountingAllocator<Z>>;
	using Entry = std::pair<const Pair, Values>;
	using Map = std::unordered_map<Pair, Values, Pair::Hash, std::equal_to<Pair>, CountingAllocator<Entry>>;
	using NativeMap = std::unordered_map<Pair, std::vector<Z>, Pair::Hash>;
	static_assert(sizeof(Values) == sizeof(std::vector<Z>));
	static_assert(sizeof(Entry) == sizeof(NativeMap::value_type));
	static_assert(sizeof(Map) == sizeof(NativeMap));
	std::cout << "Z_bytes=" << sizeof(Z) << " Pair_bytes=" << sizeof(Pair)
		<< " vector_bytes=" << sizeof(Values) << " map_value_bytes=" << sizeof(Entry)
		<< " map_object_bytes=" << sizeof(Map) << '\n';
	std::cout << "admitted,states,pairs,buckets,value_allocations,value_bytes,bucket_allocations,bucket_bytes,node_allocations,node_bytes,state_bytes_per_pair,requested_hash_bytes_per_pair,break_even_density\n";
	for (size_t states : {size_t(4), size_t(9), size_t(26), size_t(66)}) {
		for (bool pairs : {false, true}) {
			requested = {};
			Map map;
			map.reserve(admitted);
			for (size_t i=0; i<admitted; ++i)
				map.try_emplace(Pair{i/100, i%100}, states*(pairs?2:1), Z(0));
			const size_t stateBytes = states*(pairs?2:1)*sizeof(Z);
			const double bytesPerPair = double(requested[0].bytes+requested[1].bytes+requested[2].bytes)/admitted;
			std::cout << admitted << ',' << states << ',' << pairs << ',' << map.bucket_count();
			for (const auto &counter : requested) std::cout << ',' << counter.calls << ',' << counter.bytes;
			std::cout << ',' << stateBytes << ',' << std::setprecision(12) << bytesPerPair << ',' << stateBytes/bytesPerPair << '\n';
		}
	}
}

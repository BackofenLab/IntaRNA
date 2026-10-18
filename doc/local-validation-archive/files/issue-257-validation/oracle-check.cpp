#include "SeededTheoryFixtures.h"
#include <cassert>
#include <iostream>
int main() { for(const auto & f:seeded_oracle::theoryFixtures()) { auto r=seeded_oracle::enumerate(f.model); assert(std::abs(r.z-f.expected.z)<2e-12L*f.expected.z); for(auto [p,m]:f.expected.mass) assert(std::abs(r.mass[p]-m)<2e-12L*std::max(1.L,m)); std::cout<<f.name<<" "<<r.z<<"\n"; } }

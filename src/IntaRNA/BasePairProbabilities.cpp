#include "IntaRNA/BasePairProbabilities.h"
#include "IntaRNA/PartitionArithmetic.h"
#include <algorithm>
#include <stdexcept>

namespace IntaRNA {
void BasePairProbabilities::markApproximate() {
	if (state!=Status::pending) throw std::logic_error("pair probabilities: accumulator is not pending");
	approximate=true;
}
BasePairProbabilities::BasePairProbabilities(size_t n,size_t m,bool collectSeedPairs)
	: mass(n,m,0), collectSeeds(collectSeedPairs), seeds(collectSeedPairs?n:0,collectSeedPairs?m:0,0) {}
void BasePairProbabilities::addRegion(const IndexRange & t,const IndexRange & q,
		Z_type denominator,const Matrix<Z_type> & values,const Matrix<unsigned char> * seedPairs)
{
	using A=PartitionArithmetic;
	try {
		if (state!=Status::pending) throw std::logic_error("pair probabilities: accumulator is not pending");
		if (!t.isAscending() || !q.isAscending() || t.to>=mass.size1() || q.to>=mass.size2()
				|| values.size1()!=t.to-t.from+1 || values.size2()!=q.to-q.from+1)
			throw std::invalid_argument("pair probabilities: invalid region dimensions");
		if (collectSeeds && (!seedPairs || seedPairs->size1()!=values.size1() || seedPairs->size2()!=values.size2()))
			throw std::invalid_argument("pair probabilities: missing or invalid regional seed mask");
		for (const auto & r:regions)
			if (t.from<=r.first.to && r.first.from<=t.to && q.from<=r.second.to && r.second.from<=q.to)
				throw std::invalid_argument("pair probabilities: searched regions overlap");
		const Z_type merged=A::add(z,denominator);
		for (size_t i=0;i<values.size1();++i) for(size_t j=0;j<values.size2();++j) {
			A::check(values(i,j));
			if (denominator==0 && values(i,j)!=0) throw std::logic_error("pair probabilities: nonzero mass in empty region");
			A::add(mass(t.from+i,q.from+j),values(i,j));
		}
		// Allocation and all potentially failing arithmetic precede modification.
		regions.emplace_back(t,q);
		for (size_t i=0;i<values.size1();++i) for(size_t j=0;j<values.size2();++j)
		{
			mass(t.from+i,q.from+j)+=values(i,j);
			if (collectSeeds) seeds(t.from+i,q.from+j)=(*seedPairs)(i,j)!=0;
		}
		z=merged;
	} catch (...) { fail(); throw; }
}
Matrix<Z_type> BasePairProbabilities::probabilities() const {
	using A=PartitionArithmetic;
	if (state!=Status::nonempty) throw std::logic_error("pair probabilities: no successful nonempty result");
	Matrix<Z_type> p(mass.size1(),mass.size2(),0);
	// Only roundoff-sized endpoint correction is allowed. Bound scales with
	// matrix dimensions, at 512 ulps per row/column element (tested with oracles).
	const Z_type tolerance=512*std::numeric_limits<Z_type>::epsilon()
			*Z_type(std::max(size_t(1),mass.size1()+mass.size2()));
	std::vector<Z_type> columns(mass.size2(),0);
	for(size_t i=0;i<mass.size1();++i) {
		Z_type row=0;
		for(size_t j=0;j<mass.size2();++j) {
			Z_type value=A::divide(mass(i,j),z);
			if (value>1+tolerance) throw std::logic_error("pair probabilities: pair mass exceeds denominator");
			value=std::min(Z_type(1),value);
			p(i,j)=value;row=A::add(row,value);columns[j]=A::add(columns[j],value);
		}
		if(row>1+tolerance) throw std::logic_error("pair probabilities: target exclusivity violated");
	}
	for(auto c:columns) if(c>1+tolerance) throw std::logic_error("pair probabilities: query exclusivity violated");
	return p;
}
void BasePairProbabilities::finalize() {
	try {
		if(state!=Status::pending) throw std::logic_error("pair probabilities: completion requires a pending result");
		PartitionArithmetic::check(z);
		state=z==0?Status::empty:Status::nonempty;
		if(z>0) probabilities(); // validate all entries before a writer can observe success
	} catch (...) { fail();throw; }
}
}

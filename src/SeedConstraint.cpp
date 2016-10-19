
#include "SeedConstraint.h"
#include "general.h"

#include <cmath>


/////////////////////////////////////////////////////////////////////////////

SeedConstraint::SeedConstraint(
		 const size_t bp_
		, const size_t maxUnpairedOverall_
		, const size_t maxUnpaired1_
		, const size_t maxUnpaired2_
		, const E_type maxE_
		)
 :
	  bp(bp_)
	, maxUnpairedOverall(maxUnpairedOverall_)
	, maxUnpaired1(std::min(maxUnpaired1_,maxUnpairedOverall_)) // exclude too large boundaries
	, maxUnpaired2(std::min(maxUnpaired2_,maxUnpairedOverall_)) // exclude too large boundaries
	, maxE(maxE_)
{
	if (bp < 2) throw std::runtime_error("SeedHandler() : base pair number < 2 ("+toString(bp)+")");
}

/////////////////////////////////////////////////////////////////////////////

SeedConstraint::~SeedConstraint() {
}

/////////////////////////////////////////////////////////////////////////////

size_t
SeedConstraint::
getBasePairs() const {
	return bp;
}

/////////////////////////////////////////////////////////////////////////////

E_type
SeedConstraint::
getMaxE() const {
	return maxE;
}

/////////////////////////////////////////////////////////////////////////////

size_t
SeedConstraint::
getMaxUnpaired1() const {
	return maxUnpaired1;
}

/////////////////////////////////////////////////////////////////////////////

size_t
SeedConstraint::
getMaxUnpaired2() const {
	return maxUnpaired2;
}

/////////////////////////////////////////////////////////////////////////////

size_t
SeedConstraint::
getMaxUnpairedOverall() const {
	return maxUnpairedOverall;
}

/////////////////////////////////////////////////////////////////////////////

size_t
SeedConstraint::
getMaxLength1() const {
	return getBasePairs() + getMaxUnpaired1();
}

/////////////////////////////////////////////////////////////////////////////

size_t
SeedConstraint::
getMaxLength2() const {
	return getBasePairs() + getMaxUnpaired2();
}

/////////////////////////////////////////////////////////////////////////////

std::ostream&
operator<<(std::ostream& out, const SeedConstraint& c)
{
	out <<"SeedConstraint( bp="<<c.getBasePairs()
			<<", up="<<c.getMaxUnpairedOverall()
			<<", up1="<<c.getMaxUnpaired1()
			<<", up2="<<c.getMaxUnpaired2()
			<<", E="<<c.getMaxE()
			<<")";
	return out;
}

/////////////////////////////////////////////////////////////////////////////


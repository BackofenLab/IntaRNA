
#include "IntaRNA/OutputConstraint.h"

namespace IntaRNA {

/////////////////////////////////////////////////////////////////////////////

OutputConstraint::OutputConstraint(
		  const size_t reportMax
		, const ReportOverlap reportOverlap
		, const E_type maxE
		, const E_type deltaE
		, const bool bestSeedOnly
		, const bool noLP
		, const bool needZall
		)
 :
	  reportMax(reportMax)
	, reportOverlap(reportOverlap)
	, maxE(maxE)
	, deltaE(deltaE)
	, bestSeedOnly(bestSeedOnly)
	, noLP(noLP)
	, needZall(needZall)
{
	if(deltaE < (E_type)0.0) throw std::runtime_error("OutputConstraint(deltaE="+toString(deltaE)+") not >= 0.0");
}

/////////////////////////////////////////////////////////////////////////////

OutputConstraint::~OutputConstraint()
{
}

/////////////////////////////////////////////////////////////////////////////


} // namespace

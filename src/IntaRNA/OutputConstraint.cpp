
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
		, const bool noGUend
		, const bool needZall
		, const bool needBPs
		, const E_type maxED
		)
 :
	  reportMax(reportMax)
	, reportOverlap(reportOverlap)
	, maxE(maxE)
	, deltaE(deltaE)
	, bestSeedOnly(bestSeedOnly)
	, noLP(noLP)
	, noGUend(noGUend)
	, needZall(needZall)
	, needBPs(needBPs)
	, maxED(maxED)
{
	if(deltaE < (E_type)0.0) throw std::runtime_error("OutputConstraint(deltaE="+toString(deltaE)+") not >= 0.0");
}

/////////////////////////////////////////////////////////////////////////////

OutputConstraint::~OutputConstraint()
{
}

/////////////////////////////////////////////////////////////////////////////


} // namespace

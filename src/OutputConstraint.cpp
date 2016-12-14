/*
 * OutputConstraint.cpp
 *
 *  Created on: 14.12.2016
 *      Author: Mmann
 */

#include "OutputConstraint.h"


/////////////////////////////////////////////////////////////////////////////

OutputConstraint::OutputConstraint(
		  const size_t reportMax
		, const ReportOverlap reportOverlap
		, const E_type deltaE )
 :
	  reportMax(reportMax)
	, reportOverlap(reportOverlap)
	, deltaE(deltaE)
{
	assert(deltaE >= 0.0);
}

/////////////////////////////////////////////////////////////////////////////

OutputConstraint::~OutputConstraint()
{
}

/////////////////////////////////////////////////////////////////////////////


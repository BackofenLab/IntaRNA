

#include "IntaRNA/Accessibility.h"

namespace IntaRNA {

////////////////////////////////////////////////////////////////////

const E_type Accessibility::ED_UPPER_BOUND = (E_type) E_INF;

////////////////////////////////////////////////////////////////////

std::ostream&
operator<<(std::ostream& out, const Accessibility& acc)
{
	const std::string delimiter = " ";
	// print one line for each i
	out <<"\n# ED values for "<<acc.getSequence() <<":"<<"\n";
	// store current flags
	std::ios_base::fmtflags oldFlags = out.flags();
	// set number output notation to scientific format
	out <<std::scientific;
	out.precision(2);
	// empty dummy for one number
	const std::string edPlaceholder = "        ";
	for (size_t i=0; i<acc.getSequence().size(); i++) {
		for (size_t t=0; t<i; t++) {
			out <<' ' << edPlaceholder << delimiter;
		}
		// print first (without delimiter)
		out <<acc.getED(i, i);
		// print remaining with delimiter
		for (size_t j=i+1; j<std::min(acc.getMaxLength()+i,acc.getSequence().size());j++) {
			out <<delimiter <<(acc.getED(i, j)>=0?" ":"") <<acc.getED(i, j);
		}
		out <<"\n";
	}
	// reset number output notation
	out.flags( oldFlags );
	// flush output stream
	out.flush();

	// return altered stream
	return out;
}

////////////////////////////////////////////////////////////////////

void
Accessibility::
writeRNAplfold_text( std::ostream& out, const Z_type RT, const bool writeProbs ) const
{
	// store current flags
	std::ios_base::fmtflags oldFlags = out.flags();
	// set number output notation to fixed format
	out <<std::scientific;
	out.precision(6);

	// header
	if (writeProbs) {
		out <<"#unpaired probabilities"
			<<'\n';
	} else {
		out <<"#ensemble delta energy to unpair a region ED"
			<<'\n';
	}

	// length information
	out <<" #i$	l=";
	for (size_t l=1; l<=std::min(getSequence().size(),getMaxLength()); l++) {
		out <<l <<'\t';
	}
	out <<'\n';

	// write values for each possible window end j
	for (size_t j=0; j<getSequence().size(); j++) {
		// print end of window = j
		out <<(j+1) <<'\t';
		size_t maxL = std::min(j+1,getMaxLength());
		// for each (increasing) window length, print value
		for ( size_t l = 1; l <= maxL; l++ ) {
			if (writeProbs) {
				// check if ED value is not infinity
				if (E_isINF( getED(j+1-l, j) ) && !(getED(j+1-l,j) < ED_UPPER_BOUND)) {
					out <<0 <<'\t';
				} else {
					// compute unpaired probability
					double value = (double)( Z_exp( - E_2_Z(getED(j+1-l, j)) / RT ) );
					// check for nan result of conversion
					if ( value != value ) {
						out <<0 <<'\t';
					} else {
						// print unpaired probability
						out <<value <<'\t';
					}
				}
			} else {
				// write ED value in kcal/mol (ensure not printing infinity)
				out <<E_2_Ekcal(std::min<E_type>( E_MAX, getED(j+1-l, j) )) <<'\t';
			}
		}
		// print NA for remaining entries
		for ( size_t l = maxL+1; l<=getMaxLength(); l++ ) {
			out <<"NA\t";
		}
		// line break
		out <<'\n';
	}
	// flush output stream
	out.flush();
	// reset number output notation
	out.flags( oldFlags );
}

////////////////////////////////////////////////////////////////////

IndexRangeList
Accessibility::
decomposeByMaxED( const size_t maxRangeLength, const size_t winSize, const size_t minRangeLength ) const
{
#if INTARNA_IN_DEBUG_MODE
	if (winSize == 0) throw std::runtime_error("Accessibility::decomposeByMaxED() : winSize == 0");
	if (winSize >= maxRangeLength) throw std::runtime_error("Accessibility::decomposeByMaxED() : winSize >= maxRangeLength");
#endif

	// container for the resulting ranges of accessible regions
	IndexRangeList accessibleRanges;
	// insert full sequence range
	const size_t seqLastIdx = getSequence().size()-1;
	accessibleRanges.push_back( IndexRange(0,seqLastIdx) );

	for (size_t rIdx = 0; rIdx < accessibleRanges.size(); /*increment in loop*/) {
		// access to current range
		IndexRange & curRange = accessibleRanges.get(rIdx);

		// check if range is below minimal length
		if (curRange.to-curRange.from+1 < minRangeLength) {
			// remove curRange from list
			IndexRangeList::iterator curRangeIter = accessibleRanges.begin();
			std::advance( curRangeIter, rIdx );
			accessibleRanges.erase( curRangeIter );
			// skip further processing since we erased the range
			continue;
		}

		// check if to be decomposed
		if (curRange.to-curRange.from >= maxRangeLength) {
			// get boundaries for window-screening
			const size_t minIdx = curRange.from;
			const size_t maxIdx = curRange.to - winSize + 1;
			// find window with highest ED
			size_t maxEdIdx = curRange.from;
			E_type maxEd = getED(maxEdIdx,maxEdIdx+winSize-1);
			for (size_t i = minIdx; i <= maxIdx; i++ ) {
				const size_t winSizeEnd = i+winSize-1;
				// check if we found window with higher ED
				if (getED(i,winSizeEnd) > maxEd) {
					// set to middle of window
					maxEdIdx = i;
					// update maximum found so far
					maxEd = getED(i,winSizeEnd);
				}
			}
			// decompose
			if (maxEdIdx == curRange.from) {
				// prune start window
				curRange.from = std::min(maxEdIdx+winSize,curRange.to);
			} else if (maxEdIdx >= maxIdx) {
				// prune end
				curRange.to = std::max(maxEdIdx-1,curRange.from);
			} else if (maxEdIdx - curRange.from < minRangeLength) {
				// skip range in front of maxEdIdx since too short
				curRange.from = std::min(maxEdIdx+winSize,curRange.to);
			} else {
				// range after decomposition cut-point
				IndexRange newRange(std::min(maxEdIdx+winSize,curRange.to), curRange.to);
				// update current range
				curRange.to=std::max(maxEdIdx-1,curRange.from);
				// store new range if long enough
				if (newRange.to-newRange.from+1 >= minRangeLength) {
					accessibleRanges.insert(newRange);
				}
			}
		} else {
			// proceed to next range
			rIdx++;
		}
	}

	// return final ranges
	return accessibleRanges;
}


////////////////////////////////////////////////////////////////////

void
Accessibility::
decomposeByMaxED( IndexRangeList & ranges, const E_type maxED ) const
{
	// check if there is an upper bound given; if not stop working
	if (maxED >= ED_UPPER_BOUND) {
		return;
	}

	// the range list to fill
	IndexRangeList out;

	// decompose each range individually
	for (auto range = ranges.begin(); range != ranges.end(); range++) {

		size_t lastStart = range->to +1;
		for (size_t i= range->from; i <= range->to; i++) {
			if (E_isINF(getED(i,i)) || (getED(i,i) > maxED && !E_equal(getED(i,i),maxED))) {
				// check if end of range found and to be stored
				if (lastStart < i) {
					out.push_back(IndexRange(lastStart,i));
				}
				lastStart = range->to +1;
			} else {
				// check if beginning of new range
				if (lastStart > i) {
					lastStart = i;
				}
			}
		}
		// handle last open range
		if (lastStart <= range->to) {
			out.push_back(IndexRange(lastStart,range->to));
		}
	}

	// replace input ranges with final decomposed list
	ranges = out;
}

////////////////////////////////////////////////////////////////////

} // namespace

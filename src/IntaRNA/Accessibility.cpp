

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
writeRNAplfold_text( std::ostream& out, const E_type RT, const bool writeProbs ) const
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
					double value = ( std::exp( - getED(j+1-l, j) / RT ) );
					// check for nan result of conversion
					if ( value != value ) {
						out <<0 <<'\t';
					} else {
						// print unpaired probability
						out <<value <<'\t';
					}
				}
			} else {
				// write ED value (ensure not printing infinity)
				out <<std::min<E_type>( std::numeric_limits<E_type>::max(), getED(j+1-l, j) ) <<'\t';
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

void
Accessibility::
compute_accessible_ranges( IndexRangeList& ranges, const size_t max_seq_length) const
{
	if(getSequence().size() < max_seq_length){
		ranges.push_back(IndexRange(0,getSequence().size()-1));
	}
	else{
		compute_accessible_ranges_recursive( ranges, max_seq_length, ED_UPPER_BOUND, 0, getSequence().size());
	}
}

////////////////////////////////////////////////////////////////////

void
Accessibility::
compute_accessible_ranges_recursive( IndexRangeList& ranges, const size_t max_seq_length, E_type threshold, int i, int j) const
{

	const int l = 50; // Arbitrarily chosen window
	int from = -1;
	int to = -1;
	E_type ED = ED_UPPER_BOUND;

	for (i; i<j+1-l; i++) {

		if(std::min<E_type>( ED_UPPER_BOUND, getED(i, i+l) ) < threshold){
			if(from == -1 && i > to){
				from = i;
				to = i+l;
	 			ED = getED(i, i+l);
			}
			else if(from != -1 && i <= (to-l/2)){
				
				if(ED < getED(i, i+l) && (i+l-to) == 1){
					ED = getED(i, i+l);
				}
				to = i+l;
			}
			if(from != -1 && (i == j-l || (i > (to-l/2) && i <= to))){				
				if((to-from) > max_seq_length){
					compute_accessible_ranges_recursive(ranges, max_seq_length, ED, from, to);
				}
				else{
					ranges.push_back(IndexRange(from,to));
				}
				from = -1;
				ED = ED_UPPER_BOUND;
			}
		}
		else if(from != -1 && (i == j-l || i > to)){			
			if((to-from) > max_seq_length){
				compute_accessible_ranges_recursive(ranges, max_seq_length, ED, from, to);
			}
			else{
				ranges.push_back(IndexRange(from,to));
			}
			from = -1;
			ED = ED_UPPER_BOUND;
		}
	}

	return;
}

////////////////////////////////////////////////////////////////////

} // namespace

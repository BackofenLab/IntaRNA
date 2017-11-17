

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

} // namespace

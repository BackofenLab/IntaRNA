

#include "Accessibility.h"

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


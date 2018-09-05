
#include <limits>

#include "IntaRNA/IndexRange.h"

namespace IntaRNA {

////////////////////////////////////////////////////////////////

const size_t IndexRange::NA_INDEX = std::numeric_limits<size_t>::max();

////////////////////////////////////////////////////////////////

const size_t IndexRange::LAST_INDEX = IndexRange::NA_INDEX-1;

////////////////////////////////////////////////////////////////

const boost::regex IndexRange::regex("^(\\d|[123456789]\\d*)-(\\d|[123456789]\\d*)$");

////////////////////////////////////////////////////////////////

std::vector<IndexRange>
IndexRange::
overlappingWindows(const size_t windowWidth, const size_t windowOverlap) const
{
#if INTARNA_IN_DEBUG_MODE
	// minimal width check
	if (windowWidth <= 2) {
		throw std::runtime_error("IndexRange::overlappingWindows("+toString(windowWidth)+","+toString(windowOverlap)+") : "
				+"window width must be at least 2");
	}
	// ensure width exceeds overlap
	if (windowWidth <= windowOverlap) {
		throw std::runtime_error("IndexRange::overlappingWindows("+toString(windowWidth)+","+toString(windowOverlap)+") : "
				+"window width must be larger than the overlap");
	}
	// ensure ascending range
	if (isDescending()) {
		throw std::runtime_error("IndexRange::overlappingWindows() : only implemented for ascending ranges but called on "+toString(*this));
	}
#endif

	// check if range is smaller than window width
	if ((to - from + 1) <= windowWidth)
	{
		// return input index range
		return std::vector<IndexRange>(1,IndexRange(from,to));
	}

	// size_t numberOfWindows = ceil(double(to - from - windowOverlap + 1) / (windowWidth - windowOverlap));
	size_t x = to - from + 1 - windowOverlap;
	size_t y = windowWidth - windowOverlap;

#if INTARNA_IN_DEBUG_MODE
	// avoid overflows
	if (((std::numeric_limits<size_t>::max)() - x) < y)
	{
		throw std::runtime_error("IndexRange::overlappingWindows("+toString(windowWidth)+","+toString(windowOverlap)+") : "
			+ "an overflow occurred when calculating the number of windows for range "+toString(*this));
	}
#endif

	// allocate window container
	size_t numberOfWindows = (x + y - 1) / y;
	std::vector<IndexRange> windows = std::vector<IndexRange>(numberOfWindows);

	// compute windows and store
	int i = -1;
	size_t currentIndex = from;
	while (currentIndex + windowOverlap - 1 < to)
	{
		i++;
		windows[i] = IndexRange(currentIndex, currentIndex + windowWidth - 1);
		currentIndex = windows[i].to - windowOverlap + 1;
	}

	// handle last window (cut at upper boundary of this range)
	if (windows[i].to > to)
	{
		windows[i].to = to;
	}

	return windows;
}

////////////////////////////////////////////////////////////////

} // namespace

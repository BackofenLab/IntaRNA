#ifndef INTARNA_PARTITIONARITHMETIC_H_
#define INTARNA_PARTITIONARITHMETIC_H_
#include "IntaRNA/general.h"
#include <stdexcept>

namespace IntaRNA {
/** Checked nonnegative arithmetic for pair-probability partitions.
 * Structural zero is exactly zero. Positive subnormals are rejected because
 * their relative precision is not guaranteed. Ordinary sum rounding remains.
 */
class PartitionArithmetic {
public:
	/** @return x; throws on negative, NaN, infinite or positive subnormal values */
	static Z_type check(const Z_type x);
	/** @return checked sum of two nonnegative values */
	static Z_type add(const Z_type a, const Z_type b);
	/** @return checked product; positive underflow is an error */
	static Z_type multiply(const Z_type a, const Z_type b);
	/** @return checked quotient; denominator must be positive */
	static Z_type divide(const Z_type a, const Z_type b);
	/** @return checked exponential; a finite exponent must yield a positive normal value */
	static Z_type exp(const Z_type exponent);
};
inline Z_type PartitionArithmetic::check(const Z_type x) {
	if (!(x>=0 && x<=std::numeric_limits<Z_type>::max())
			|| (x>0 && x<std::numeric_limits<Z_type>::min()))
		throw std::range_error("partition arithmetic: numerical range exceeded");
	return x;
}
inline Z_type PartitionArithmetic::add(const Z_type a,const Z_type b) {
	return check(check(a)+check(b));
}
inline Z_type PartitionArithmetic::multiply(const Z_type a,const Z_type b) {
	check(a); check(b);
	if (a==0 || b==0) return 0;
	const Z_type r=check(a*b);
	if (r==0) throw std::range_error("partition arithmetic: product underflow");
	return r;
}
inline Z_type PartitionArithmetic::divide(const Z_type a,const Z_type b) {
	check(a); check(b);
	if (b==0) throw std::range_error("partition arithmetic: zero denominator");
	const Z_type r=check(a/b);
	if (a>0 && r==0) throw std::range_error("partition arithmetic: quotient underflow");
	return r;
}
inline Z_type PartitionArithmetic::exp(const Z_type exponent) {
	const Z_type r=check(Z_exp(exponent));
	if (r==0) throw std::range_error("partition arithmetic: exponential underflow");
	return r;
}
}
#endif

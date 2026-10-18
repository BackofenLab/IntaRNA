#include <IntaRNA/Matrix.h>
static_assert(INTARNA_USE_STD_MDSPAN == 0);
int main() {
    IntaRNA::Matrix<int> dense(2,3); dense(1,2)=42;
    IntaRNA::UpperTriangularMatrix<int> triangle(3,3); triangle(1,2)=7;
    IntaRNA::UpperBandedMatrix<int> band(3,3,0,1); band(1,2)=9;
    return dense(1,2)==42 && triangle(1,2)==7 && band(1,2)==9 ? 0 : 1;
}

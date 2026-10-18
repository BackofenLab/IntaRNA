#include <IntaRNA/Matrix.h>
#include <IntaRNA/RnaSequence.h>
#include <utility>
static_assert(INTARNA_USE_STD_MDSPAN == 0);
int main() {
    IntaRNA::RnaSequence rna("test", "ACGU");
    IntaRNA::Matrix<int> dense(2,3); dense(1,2)=42;
    IntaRNA::UpperTriangularMatrix<int> triangle(4,4); triangle(2,3)=7;
    IntaRNA::UpperBandedMatrix<int> band(4,4,0,1); band(2,3)=9;
    return rna.size()==4 && dense(1,2)==42 && triangle(2,3)==7
        && band(2,3)==9 && std::as_const(triangle)(3,2)==0
        && std::as_const(band)(1,3)==0 ? 0 : 1;
}

#include "IntaRNA/general.h"
#include <iomanip>
#include <iostream>
#include <sstream>
int main() {
    using namespace IntaRNA;
    for (const char* text : {"4.760000e+00", "4.890000e+00", "4.770000e+00", "9.800000e+00", "1.086000e+01"}) {
        double value = 0; // AccessibilityFromStream reads into double curVal.
        std::istringstream(text) >> value;
        std::cout << text << '\t' << std::setprecision(17) << value*100 << '\t' << Ekcal_2_E(value) << '\n';
    }
}

#include <IntaRNA/AccessibilityVrna.h>
#include <IntaRNA/NussinovHandler.h>
#include <IntaRNA/PredictorMfeEns2dSeedExtension.h>
int main() {
    IntaRNA::Matrix<int> dense(2, 3, 7);
    IntaRNA::UpperBandedMatrix<int> band(3, 3, 0, 1);
    band(1, 2) = dense(1, 2);
    IntaRNA::NussinovHandler::Z2dMatrix triangle(3, 3);
    triangle(0, 2) = band(1, 2);
    return triangle(0, 2) == 7 ? 0 : 1;
}

#include <IntaRNA/VrnaHandler.h>
#include <IntaRNA/RnaSequence.h>
#include <IntaRNA/AccessibilityVrna.h>
INITIALIZE_EASYLOGGINGPP
int main() {
    IntaRNA::VrnaHandler handler;
    IntaRNA::RnaSequence sequence("test", "GGGGAAAACCCC");
    IntaRNA::AccessibilityVrna accessibility(sequence, 6, nullptr, handler, 10);
    return accessibility.getED(0, 2) < IntaRNA::E_INF ? 0 : 1;
}

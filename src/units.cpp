#include "units.h"

double Faunus::PhysicalConstants::temperature = 298.15;

std::string Faunus::u8::bracket(const std::string &s) {
    return "\u27e8" + s + "\u27e9";
}


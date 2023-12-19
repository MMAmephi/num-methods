#include <iostream>
#include <limits>

int main() {
    /* 1. Epsilon */

    std::cout << "Epsilon:\n";
    // float32
    float eps32 = 1.0;
    while ((float)1.0 + eps32 > (float)1.0) {
        eps32 /= (float)2.0;
    }
    std::cout << static_cast<float>(eps32 * 2) << std::endl;

    // double64
    double eps64 = 1.0;
    while (1.0 + eps64 > 1.0) {
        eps64 /= 2.0;
    }
    std::cout << static_cast<double>(eps64 * 2) << std::endl;

    // long double128
    long double eps128 = 1.0;
    while (1.0 + eps128 > 1.0) {
        eps128 /= 2.0;
    }
    std::cout << eps128 * 2.0 << std::endl;

    /* 2. Dwarf */

    std::cout << "Dwarf:\n";
    // float32
    float dwarf32 = 1.0;
    float temp32 = 0.0;
    while (dwarf32 > (float)0.0) {
        temp32 = dwarf32;
        dwarf32 /= (float)2.0;
    }
    std::cout << temp32 << std::endl;

    // double 128
    double dwarf64 = 1.0;
    double temp64 = 0.0;
    while (dwarf64 > 0.0) {
        temp64 = dwarf64;
        dwarf64 /= 2.0;
    }
    std::cout << temp64 << std::endl;

    // long double 128
    long double dwarf128 = 1.0;
    long double temp128 = 0.0;
    while (dwarf128 > 0.0) {
        temp128 = dwarf128;
        dwarf128 /= 2.0;
    }
    std::cout << temp128 << std::endl;

    /* 3. Max */

    std::cout << "Max number:\n";

    // float32
    float max32 = 1.0;
    float tempMax32 = 0.0;
    while (max32 != std::numeric_limits<float>::infinity()) {
        tempMax32 = max32;
        max32 *= 2.0;
    }
    std::cout << tempMax32 << std::endl;

    // double 64
    double max64 = 1.0;
    double tempMax64 = 0.0;
    while (max64 != std::numeric_limits<double>::infinity()) {
        tempMax64 = max64;
        max64 *= 2.0;
    }
    std::cout << tempMax64 << std::endl;

    // long double 128
    long double max128 = 1.0;
    long double tempMax128 = 0.0;
    while (max128 != std::numeric_limits<long double>::infinity()) {
        tempMax128 = max128;
        max128 *= 2.0;
    }
    std::cout << tempMax128 << std::endl;

    return 0;
}
#pragma once

#include <vector>
#include <iostream>

class VectorPrinter {
public:

    static void printVector(const std::vector<double>& v) {
        std::cout << "[" << v[0];
        for (size_t i = 1; i < v.size(); i++)
        {
            std::cout << ", " << v[i];
        }
        std::cout << "]\n";
    }

private:
    VectorPrinter() {}
    ~VectorPrinter() {}

};

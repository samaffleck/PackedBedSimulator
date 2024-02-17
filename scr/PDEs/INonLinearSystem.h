#pragma once

#include <vector>
#include <iostream>

#include "../VectorFunctions/VectorNorm.h"


class INonLinearSystem {
public:
    INonLinearSystem(
        std::vector<double> _x,
        int _order) :
        x(_x),
        order(_order)
    {
        e.resize(x.size());
        for (int i = 0; i < order; i++)
        {
            xPrev.push_back(_x);            // Copy initial vector to the previous vector for all orders.
        }
    }
    ~INonLinearSystem() {}

    int order{};
    int itterations{};
    double error{};

    std::vector<double> x;                     // Vector of scalar elements
    std::vector<double> e;                     // Function result
    std::vector<std::vector<double>> xPrev;    // Vectors at previous time steps depending on the order...

    virtual void gaussSeidel(const double& dt) = 0;
    virtual double evaluateError(const double& dt) = 0;

    void innerItteration(const int& maxItterations, const double& tolerance, const double& timeStep) {
        error = tolerance * 100;
        itterations = 0;

        while (error > tolerance && itterations < maxItterations) {
            gaussSeidel(timeStep);
            error = evaluateError(timeStep);
            itterations++;
        }

        std::cout << "error = " << error << "\tInner Itterations = " << itterations << "\n";

    }

};

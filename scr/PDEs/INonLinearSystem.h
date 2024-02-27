#pragma once

#include <vector>
#include <iostream>

#include "../VectorFunctions/VectorNorm.h"
#include "../BoundaryConditions/IBoundaryCondition.h"
#include "../Steps/IStep.h"

class INonLinearSystem {
public:

    INonLinearSystem() {}
    INonLinearSystem(
        std::vector<double> _x,
        int _order,
        double _length,
        int _numberOfcells) :
        length(_length),
        numberOfcells(_numberOfcells),
        x(_x),
        order(_order)
    {
        dx = length / _numberOfcells;
        e.resize(x.size());
        for (int i = 0; i < order; i++)
        {
            xPrev.push_back(_x);            // Copy initial vector to the previous vector for all orders.
        }
        rhs.resize(x.size());
    }
    ~INonLinearSystem() {}


    int order{};
    int itterations{};
    double error{};

    double length{};
    int numberOfcells{};
    double dx{};
    const double R = 8.314;

    std::vector<double> x;                     // Vector of scalar elements
    std::vector<double> e;                     // Function result
    std::vector<std::vector<double>> xPrev;    // Vectors at previous time steps depending on the order...
    std::vector<double> rhs{};

    IStep* step = nullptr;

    virtual void updateRHS(const double& dt) = 0;

    void gaussSeidel(const double& dt);
    double evaluateError(const double& dt);
    void innerItteration(const int& maxItterations, const double& tolerance, const double& timeStep);

};

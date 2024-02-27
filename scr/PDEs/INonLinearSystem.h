#pragma once

#include <vector>
#include <iostream>

#include "../VectorFunctions/VectorNorm.h"
#include "../BoundaryConditions/IBoundaryCondition.h"
#include "../Steps/IStep.h"
#include "../SystemObjects/PackedBed.h"

class PackedBed;

class INonLinearSystem {
public:

    INonLinearSystem(
        PackedBed* _bed,
        std::vector<double>& _x,
        int _order) :
        bed(_bed),
        x(_x),
        order(_order)
    {
        e.resize(bed->numberOfCells);
        for (int i = 0; i < order; i++)
        {
            xPrev.push_back(_x);            // Copy initial vector to the previous vector for all orders.
        }
        rhs.resize(bed->numberOfCells);
    }
    ~INonLinearSystem() {}


    int order{};
    int itterations{};
    double error{};

    const double R = 8.314;

    std::vector<double>& x;                    // Vector of scalar elements
    std::vector<double> e;                     // Function result
    std::vector<std::vector<double>> xPrev;    // Vectors at previous time steps depending on the order...
    std::vector<double> rhs{};

    PackedBed* bed = nullptr;

    virtual void updateRHS(const double& dt) = 0;

    void gaussSeidel(const double& dt);
    double evaluateError(const double& dt);
    void innerItteration(const int& maxItterations, const double& tolerance, const double& timeStep);

};

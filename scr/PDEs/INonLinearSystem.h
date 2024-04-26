#pragma once

#include <vector>
#include <iostream>

#include "../VectorFunctions/VectorNorm.h"
#include "../BoundaryConditions/IBoundaryCondition.h"
#include "../Steps/IStep.h"
#include "../SystemObjects/PackedBed.h"
#include "../Solver/ThomasAlgoritms.h"

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
        e.resize(_x.size());
        x_LastItter.resize(_x.size());
        x_step.resize(_x.size());
        Ae.resize(_x.size());
        Aw.resize(_x.size());
        Ao.resize(_x.size());
        So.resize(_x.size());
        for (int i = 0; i < order; i++)
        {
            xPrev.push_back(_x);            // Copy initial vector to the previous vector for all orders.
        }
        x_LastItter = _x;
        rhs.resize(_x.size());
        TDMA.initialise((int)_x.size());
        dampingFactor = 1e-2 * bed->dx * bed->dx;
    }
    ~INonLinearSystem() {}


    int order{};
    int itterations{};
    double error{};

    const double R = 8.314;

    std::vector<double>& x;                     // Vector of scalar elements
    std::vector<double> x_LastItter;            // Vector of scalar elements
    std::vector<double> x_step;                 // Vector of steps
    std::vector<double> e;                      // Function result
    std::vector<std::vector<double>> xPrev;     // Vectors at previous time steps depending on the order...
    std::vector<double> rhs{};

    std::vector<double> Ae;                     // East node link coefficient
    std::vector<double> Aw;                     // West node link coefficient
    std::vector<double> Ao;                     // Central node link coefficient
    std::vector<double> So;                     // Central node source coefficient

    double dampingFactor;                 // Damping factor

    PackedBed* bed = nullptr;

    ThomasAlgorithm TDMA;

    virtual void updateRHS(const double& dt) = 0;
    virtual void updateLinkCoefficients(const double& dt) = 0;

    void gaussSeidel(const double& alpha);
    double evaluateError();
    void innerItteration(const int& maxItterations, const double& tolerance, const double& timeStep);

};

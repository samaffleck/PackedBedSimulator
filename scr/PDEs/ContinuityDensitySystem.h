#pragma once

#include "INonLinearSystem.h"

class ContinuityDensitySystem : public INonLinearSystem {
public:

    ContinuityDensitySystem() {}
    ContinuityDensitySystem(
        double _kappa,
        double _visocisty,
        double _et,
        std::vector<double> _x,
        int _order, 
        int _numberOfCells,
        double _length) :
        kappa(_kappa),
        vis(_visocisty),
        et(_et),
        INonLinearSystem(_x, _order, _length, _numberOfCells) 
    {}
    ~ContinuityDensitySystem() {}
    ContinuityDensitySystem& operator=(const ContinuityDensitySystem& other)
    {
        if (this != &other) {
            length = other.length;
            dx = other.dx;
            numberOfcells = other.numberOfcells;
            order = other.order;
            itterations = other.itterations;
            error = other.error;
            x = other.x;
            e = other.e;
            xPrev = other.xPrev;
            rhs = other.rhs;
            step = other.step;

            kappa = other.kappa;
            vis = other.vis;
            et = other.et;
        }
        return *this;
    }

    double kappa{};
    double vis{};
    double et{};        // Total voidage

    INonLinearSystem* uSystem = nullptr;
    INonLinearSystem* tSystem = nullptr;

    void updateRHS(const double& dt) override;

};

#pragma once

#include "INonLinearSystem.h"


class ContinuityVelocitySystem : public INonLinearSystem {
public:

    ContinuityVelocitySystem() {}
    ContinuityVelocitySystem(
        double _kappa,
        double _viscosity,
        std::vector<double> _x,
        int _order,
        int _numberOfCells,
        double _length) :
        kappa(_kappa),
        vis(_viscosity),
        INonLinearSystem(_x, _order, _length, _numberOfCells) {}
    ~ContinuityVelocitySystem() {}
    ContinuityVelocitySystem& operator=(const ContinuityVelocitySystem& other)
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
        }
        return *this;
    }


    double kappa{};
    double vis{};

    INonLinearSystem* cSystem = nullptr;
    INonLinearSystem* tSystem = nullptr;

    void updateRHS(const double& dt) override;

};

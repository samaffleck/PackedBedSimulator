#pragma once

#include "INonLinearSystem.h"


class ContinuityVelocitySystem : public INonLinearSystem {
public:

    ContinuityVelocitySystem(
        std::vector<double> _x,
        int _order,
        int _numberOfCells,
        double _length) :
        INonLinearSystem(_x, _order, _length, _numberOfCells) {}
    ~ContinuityVelocitySystem() {}

    double kappa = 9e-9;
    double vis = 2e-5;
    double const1 = kappa / (vis * dx);

    INonLinearSystem* pSystem = nullptr;

    void updateRHS(const double& dt) override;

};

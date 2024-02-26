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
        INonLinearSystem(_x, _order, _length, _numberOfCells) 
    {
        const1 = kappa / (vis * dx);
    }
    ~ContinuityVelocitySystem() {}

    double kappa{};
    double vis{};
    double const1{};

    INonLinearSystem* pSystem = nullptr;

    void updateRHS(const double& dt) override;

};

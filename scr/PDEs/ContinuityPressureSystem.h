#pragma once

#include "INonLinearSystem.h"

class ContinuityPressureSystem : public INonLinearSystem {
public:

    ContinuityPressureSystem() {}
    ContinuityPressureSystem(
        double _kappa,
        double _visocisty,
        std::vector<double> _x,
        int _order, 
        int _numberOfCells,
        double _length) :
        kappa(_kappa),
        vis(_visocisty),
        INonLinearSystem(_x, _order, _length, _numberOfCells) 
    {
        const1 = kappa / (3 * vis * dx * dx);
    }
    ~ContinuityPressureSystem() {}

    double kappa{};
    double vis{};
    double const1{};

    INonLinearSystem* uSystem = nullptr;

    void updateRHS(const double& dt) override;

};

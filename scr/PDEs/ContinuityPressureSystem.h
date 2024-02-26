#pragma once

#include "INonLinearSystem.h"

class ContinuityPressureSystem : public INonLinearSystem {
public:

    ContinuityPressureSystem(
        std::vector<double> _x,
        int _order, 
        int _numberOfCells,
        double _length) :
        INonLinearSystem(_x, _order, _length, _numberOfCells) {}
    ~ContinuityPressureSystem() {}

    double kappa = 9e-9;
    double vis = 2e-5;
    double const1 = kappa / (3 * vis * dx * dx);

    INonLinearSystem* uSystem = nullptr;

    void updateRHS(const double& dt) override;

};

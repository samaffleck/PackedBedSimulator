#pragma once

#include "INonLinearSystem.h"


class DiffusionSystem : public INonLinearSystem {
public:

    DiffusionSystem(
        std::vector<double> _x,
        int _order) :
        INonLinearSystem(_x, _order){}
    ~DiffusionSystem() {}

    double dx = 0.002;
    double dxSquared = dx * dx;
    double D = 0.5;
    double D_dxSquaredInv = D / dxSquared;
    double qSource = 1e6;

    void updateRHS(const double& dt) override;

};

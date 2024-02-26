#pragma once

#include "INonLinearSystem.h"


class DiffusionSystem : public INonLinearSystem {
public:

    DiffusionSystem(
        std::vector<double> _x,
        int _order,
        int _numberOfCells,
        double _length) :
        INonLinearSystem(_x, _order, _length, _numberOfCells) {}
    ~DiffusionSystem() {}

    double dxSquared = dx * dx;
    double D = 0.5;
    double D_dxSquaredInv = D / dxSquared;
    double qSource = 1e6;

    void updateRHS(const double& dt) override;

};

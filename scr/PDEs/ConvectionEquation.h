#pragma once

#include "INonLinearSystem.h"

class ConvectionEquation : public INonLinearSystem {
public:
    ConvectionEquation(std::vector<double> _x,
        int _order,
        int _numberOfCells,
        double _length) :
        INonLinearSystem(_x, _order, _length, _numberOfCells) {}
    ~ConvectionEquation() {}

    double dxInv = 1 / dx;
    double u = 35.27;
    double* time = nullptr;

    void updateRHS(const double& dt) override;

};

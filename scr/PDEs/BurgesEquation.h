#pragma once

#include "INonLinearSystem.h"

class BurgesEquation : public INonLinearSystem {

public:
	BurgesEquation(std::vector<double> _x,
        int _order,
        int _numberOfCells,
        double _length) :
        INonLinearSystem(_x, _order, _length, _numberOfCells) {}
    ~BurgesEquation(){}

    double dxInv = 1 / dx;

    void updateRHS(const double& dt) override;

};

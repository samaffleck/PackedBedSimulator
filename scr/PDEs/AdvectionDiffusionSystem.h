#pragma once

#include "INonLinearSystem.h"


class AdvectionDiffusionSystem : public INonLinearSystem {
public:

    AdvectionDiffusionSystem() {}
    AdvectionDiffusionSystem(
        std::vector<double> _x,
        int _order,
        int _numberOfCells,
        double _length) :
        INonLinearSystem(_x, _order, _length, _numberOfCells) {}
    ~AdvectionDiffusionSystem() {}
    AdvectionDiffusionSystem& operator=(const AdvectionDiffusionSystem& other)
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

            dxSquared = other.dxSquared;
            D = other.D;
            D_dxSquaredInv = other.D_dxSquaredInv;
            qSource = other.qSource;
        }
        return *this;
    }


    double dxSquared = dx * dx;
    double D = 0.5;
    double D_dxSquaredInv = D / dxSquared;
    double qSource = 1e6;

    void updateRHS(const double& dt) override;

};

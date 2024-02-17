#pragma once

#include "INonLinearSystem.h"


class ContinuityVelocitySystem : public INonLinearSystem {
public:

    ContinuityVelocitySystem(
        std::vector<double> _x,
        int _order) :
        INonLinearSystem(_x, _order) {}
    ~ContinuityVelocitySystem() {}

    double dx = 0.05;
    double kappa = 9e-9;
    double vis = 2e-5;

    std::vector<double>* p = nullptr; // Connecting variable
    std::vector<double>* pPrev = nullptr; // Connecting variable

    void gaussSeidel(const double& dt) override {

        x[0] = -(kappa / (vis * dx)) * ((*p)[1] - (*p)[0]);

        for (size_t i = 1; i < x.size() - 1; i++)
        {
            x[i] = -(kappa / (vis * dx * 2)) * ((*p)[i + 1] - (*p)[i - 1]);
        }

        x[x.size() - 1] = -(kappa / (vis * dx)) * ((*p)[x.size() - 1] - (*p)[x.size() - 2]);

    }

    double evaluateError(const double& dt) override {

        error = 0;

        e[0] = -x[0] - (kappa / (vis * dx)) * ((*p)[1] - (*p)[0]);

        for (size_t i = 1; i < x.size() - 1; i++)
        {
            e[i] = -x[i] - (kappa / (vis * dx * 2)) * ((*p)[i + 1] - (*p)[i - 1]);
        }

        e[e.size() - 1] = -x[x.size() - 1] - (kappa / (vis * dx)) * ((*p)[x.size() - 1] - (*p)[x.size() - 2]);

        error = VectorNorm::L2Norm(e);
        return error;
    }

};

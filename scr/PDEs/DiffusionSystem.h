#pragma once

#include "INonLinearSystem.h"


class DiffusionSystem : public INonLinearSystem {
public:

    DiffusionSystem(
        std::vector<double> _x,
        int _order) :
        INonLinearSystem(_x, _order)
    {
        rhs.resize(x.size());
    }
    ~DiffusionSystem() {}

    double dx = 0.05;
    double dxSquared = dx * dx;
    double D = 1e-3;
    std::vector<double> rhs{};

    void gaussSeidel(const double& dt) override {

        updateRHS(dt);
        x = rhs;

    }

    double evaluateError(const double& dt) override {

        error = 0;
        updateRHS(dt);
        for (size_t i = 0; i < e.size(); i++)
        {
            e[i] = rhs[i] - x[i];
        }
        error = VectorNorm::L2Norm(e);
        return error;

    }

    void updateRHS(const double& dt) {

        rhs[0] = (
            6 * dxSquared * xPrev[0][0] +
            16 * D * dt * 300 +
            x[1] * (8 * dt * D)
            ) / (
                6 * dxSquared +
                24 * dt * D
                );

        for (size_t i = 1; i < x.size() - 1; i++)
        {
            rhs[i] = (
                2 * dxSquared * xPrev[0][i] +
                2 * D * dt * x[i + 1] +
                2 * D * dt * x[i - 1]
                ) / (
                    2 * dxSquared +
                    4 * D * dt
                    );
        }

        rhs[rhs.size() - 1] = (
            2 * dxSquared * xPrev[0][x.size() - 1] +
            2 * D * dt * x[x.size() - 2]
            ) / (
                2 * dxSquared +
                2 * D * dt
                );

    }

};

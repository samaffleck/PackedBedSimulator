#include "INonLinearSystem.h"


void INonLinearSystem::gaussSeidel(const double& dt) {
    updateRHS(dt);
    x = rhs;
}

double INonLinearSystem::evaluateError(const double& dt) {
    error = 0;
    updateRHS(dt);
    for (size_t i = 0; i < e.size(); i++)
    {
        e[i] = rhs[i] - x[i];
    }
    error = VectorNorm::L2Norm(e);
    return error;
}


void INonLinearSystem::innerItteration(const int& maxItterations, const double& tolerance, const double& timeStep) {
    error = tolerance * 100;
    itterations = 0;

    while (error > tolerance && itterations < maxItterations) {
        gaussSeidel(timeStep);
        error = evaluateError(timeStep);
        itterations++;
    }
}

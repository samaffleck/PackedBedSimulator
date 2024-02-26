#include "BurgesEquation.h"

void BurgesEquation::updateRHS(const double& dt) {

    rhs[0] = xPrev[0][0] - (x[0] * dt / (2 * dx)) * (x[0] - inletBoundaryCondition->boundaryValue);

    for (size_t i = 1; i < x.size() - 1; i++)
    {
        rhs[i] = xPrev[0][i] - (x[i] * dt / dx) * (x[i] - x[i - 1]);
    }

    int j = rhs.size() - 1;
    rhs[j] = xPrev[0][j] - (x[j] * dt / dx) * (x[j] - x[j - 1]);

}

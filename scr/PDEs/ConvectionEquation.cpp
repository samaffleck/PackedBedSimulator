#include "ConvectionEquation.h"


void ConvectionEquation::updateRHS(const double& dt) {

    rhs[0] = xPrev[0][0] - (u * dt * 2 / (dx)) * (x[0] - inletBoundaryCondition->boundaryValue);

    for (size_t i = 1; i < x.size(); i++)
    {
        rhs[i] = xPrev[0][i] - (dt * u / dx) * (x[i] - x[i - 1]);
    }

}

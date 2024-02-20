#include "ContinuityPressureSystem.h"


void ContinuityPressureSystem::updateRHS(const double& dt) {

    //rhs[0] = (dx * xPrev[0][0] + 
    //    dt * (2 * uSystem->inletBoundaryCondition->boundaryValue - uSystem->x[0]) * (x[0] + vis * dx * uSystem->inletBoundaryCondition->boundaryValue / kappa)) / 
    //    (dx + dt * uSystem->x[0]);
    rhs[0] = xPrev[0][0] + dt * const1 * x[0] * (4 * x[1] - 12 * x[0] + 8 * inletBoundaryCondition->boundaryValue);

    for (size_t i = 1; i < x.size() - 1; i++)
    {
        //rhs[i] = (dx * xPrev[0][i] + dt * uSystem->x[i - 1] * x[i - 1]) / (dx + dt * uSystem->x[i]);
        rhs[i] = xPrev[0][i] + dt * 3 * const1 * x[i] * (x[i + 1] - 2 * x[i] + x[i - 1]);
    
    }

    //rhs[rhs.size() - 1] = xPrev[0][xPrev.size() - 1] + dt * const1 * x[x.size() - 1] * (8 * outletBoundaryCondition->boundaryValue - 12 * x[x.size() - 1] + 4 * x[x.size() - 2]);
    rhs[rhs.size() - 1] = outletBoundaryCondition->boundaryValue;

}

#include "ContinuityPressureSystem.h"


void ContinuityPressureSystem::updateRHS(const double& dt) {

    // Feed Boundary Condition
    /*
    rhs[0] = xPrev[0][0] - 
        (2 * dt / (dx)) * 
        ((x[0] * (*uSystem).x[0]) - 
        (uSystem->inletBoundaryCondition->boundaryValue * x[0] + 
            dx * vis * uSystem->inletBoundaryCondition->boundaryValue * uSystem->inletBoundaryCondition->boundaryValue / (2 * kappa)));
    */
    
    // Pressurise Boundary Condition
    rhs[0] = xPrev[0][0] -
        (dt / dx) * (0.5 * (x[0] * uSystem->x[0] + x[1] * uSystem->x[1]) -
            inletBoundaryCondition->boundaryValue * ((2 * kappa) / (dx * vis)) * (inletBoundaryCondition->boundaryValue - x[0]));
    
    for (size_t i = 1; i < x.size(); i++)
    {
        rhs[i] = xPrev[0][i] - (dt / dx) * ((x[i] * uSystem->x[i]) - (x[i - 1] * uSystem->x[i - 1]));
    }

}

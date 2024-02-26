#include "ContinuityVelocitySystem.h"

void ContinuityVelocitySystem::updateRHS(const double& dt) {

    // Feed Boundary Condition
    /*
    rhs[0] = inletBoundaryCondition->boundaryValue / 2 -
        (kappa / (2 * dx * vis)) * (pSystem->x[1] - pSystem->x[0]);
    */
    
    // Pressurise Boundary Condition
    rhs[0] = (-kappa / (dx * vis)) * (0.5 * (pSystem->x[0] + pSystem->x[1]) - 
        pSystem->inletBoundaryCondition->boundaryValue);
        
    
    for (size_t i = 1; i < x.size() - 1; i++)
    {
        rhs[i] = (-kappa / (vis * 2 * dx)) * (pSystem->x[i + 1] - pSystem->x[i - 1]);
    }


    // Feed Boundary Condition
    /*
    rhs[rhs.size() - 1] = (-kappa / (dx * vis)) * (pSystem->outletBoundaryCondition->boundaryValue -
        0.5 * (pSystem->x[rhs.size() - 1] + pSystem->x[rhs.size() - 2]));
    */
    
    // Pressurise Boundary Condition
    rhs[rhs.size() - 1] = (-kappa / (2 * dx * vis)) * (pSystem->x[rhs.size() - 1] - pSystem->x[rhs.size() - 2]);

}

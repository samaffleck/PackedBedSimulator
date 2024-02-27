#include "AdvectionDiffusionSystem.h"


void AdvectionDiffusionSystem::updateRHS(const double& dt) {
    /*
    rhs[0] = (
        6 * dxSquared * xPrev[0][0] +
        16 * D * dt * inletBoundaryCondition->boundaryValue +
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
    */


    /*
    rhs[0] = (x[1] * D_dxSquaredInv + 
        inletBoundaryCondition->boundaryValue * 2 * D_dxSquaredInv + 
        qSource + 
        xPrev[0][0] * (1 / dt)) /
        ((1 / dt) + 3 * D_dxSquaredInv);
    */

    
    
    /*->>
    for (size_t i = 1; i < x.size() - 1; i++)
    {
        rhs[i] = (x[i + 1] * D_dxSquaredInv + 
            x[i - 1] * D_dxSquaredInv + 
            qSource + 
            xPrev[0][i] * (1 / dt)) / 
            ((1 / dt) + 2 * D_dxSquaredInv);
    }
    */
    
    
    
    /*
    int j = rhs.size() - 1;
    rhs[j] = (outletBoundaryCondition->boundaryValue * 2 * D_dxSquaredInv +
        x[j - 1] * D_dxSquaredInv +
        qSource +
        xPrev[0][j] * (1 / dt)) /
        ((1 / dt) + 3 * D_dxSquaredInv);
    */

}

#include "ContinuityVelocitySystem.h"

void ContinuityVelocitySystem::updateRHS(const double& dt) {

    rhs[0] = bed->step->inletVelocityRHS(bed);
    
    for (size_t i = 1; i < x.size() - 1; i++)
    {
        rhs[i] = (-bed->kappa * R / (bed->viscosity * 2 * bed->dx)) * 
            (bed->C[i + 1] * bed->T[i + 1] - bed->C[i - 1] * bed->T[i - 1]);
    }

    rhs[rhs.size() - 1] = bed->step->outletVelocityRHS(bed);

}

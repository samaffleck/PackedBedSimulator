#include "AdvectionDiffusionSystem.h"


void AdvectionDiffusionSystem::updateRHS(const double& dt) {
    
    rhs[0] = bed->step->inletTemperatureRHS(bed, x, xPrev, dt);
    
    for (size_t i = 1; i < x.size() - 1; i++)
    {
        rhs[i] = xPrev[0][i] - 
            (dt / bed->dx) * bed->U[i] * (x[i] - x[i - 1]) +
            (dt / (bed->dx * bed->dx)) * bed->lambda * (x[i + 1] - 2 * x[i] + x[i - 1]);
    }

    rhs[rhs.size() - 1] = bed->step->outletTemperatureRHS(bed, x, xPrev, dt);

}

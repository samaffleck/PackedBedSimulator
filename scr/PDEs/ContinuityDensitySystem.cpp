#include "ContinuityDensitySystem.h"


void ContinuityDensitySystem::updateRHS(const double& dt) {

    rhs[0] = bed->step->inletDensityRHS(bed, x, xPrev, dt);

    for (size_t i = 1; i < x.size(); i++)
    {
        rhs[i] = xPrev[0][i] - (dt / (bed->et * bed->dx)) * 
            ((x[i] * bed->U[i]) - (x[i - 1] * bed->U[i - 1]));
    }

}

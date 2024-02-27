#include "ContinuityVelocitySystem.h"

void ContinuityVelocitySystem::updateRHS(const double& dt) {

    rhs[0] = step->inletVelocityRHS(kappa, dx, vis, cSystem->x, tSystem->x);
    
    for (size_t i = 1; i < x.size() - 1; i++)
    {
        rhs[i] = (-kappa * R / (vis * 2 * dx)) * (cSystem->x[i + 1] * tSystem->x[i + 1] - cSystem->x[i - 1] * tSystem->x[i - 1]);
    }

    rhs[rhs.size() - 1] = step->outletVelocityRHS(kappa, dx, vis, cSystem->x, tSystem->x, (rhs.size() - 1));

}

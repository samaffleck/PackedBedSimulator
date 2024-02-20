#include "ContinuityVelocitySystem.h"

void ContinuityVelocitySystem::updateRHS(const double& dt) {

    rhs[0] = (-const1 / 2) * (pSystem->x[1] - (pSystem->x[1] + vis * dx * inletBoundaryCondition->boundaryValue / kappa));

    for (size_t i = 1; i < x.size() - 1; i++)
    {
        rhs[i] = (-const1 / 2) * (pSystem->x[i + 1] - pSystem->x[i - 1]);
    }

    rhs[rhs.size() - 1] = (-const1 / 2) * ((2 * pSystem->outletBoundaryCondition->boundaryValue - pSystem->x[rhs.size() - 1]) - pSystem->x[rhs.size() - 2]);

}

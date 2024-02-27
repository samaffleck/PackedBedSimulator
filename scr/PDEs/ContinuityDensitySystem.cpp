#include "ContinuityDensitySystem.h"


void ContinuityDensitySystem::updateRHS(const double& dt) {

    rhs[0] = step->inletPressureRHS(xPrev, dt, dx, x, uSystem->x, tSystem->x, vis, kappa);

    for (size_t i = 1; i < x.size(); i++)
    {
        rhs[i] = xPrev[0][i] - (dt / (et * dx)) * ((x[i] * uSystem->x[i]) - (x[i - 1] * uSystem->x[i - 1]));
    }

}

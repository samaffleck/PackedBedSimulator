#include "ContinuityVelocitySystem.h"

void ContinuityVelocitySystem::updateRHS(const double& dt) {

    double d = bed->kappa / (bed->dx * bed->viscosity);

    rhs[0] = bed->step->inletVelocityRHS(bed);
    
    for (size_t i = 1; i < x.size() - 1; i++)
    {
        
        rhs[i] = d * (bed->P[i - 1] - bed->P[i]);

    }

    rhs[rhs.size() - 1] = bed->step->outletVelocityRHS(bed);

}


void ContinuityVelocitySystem::updateVelocity() {

    x_LastItter = x;

    double d = - bed->kappa / (bed->dx * bed->viscosity);
    const auto N = (int)x.size();

    x[0] = bed->step->inletVelocityRHS(bed);

    for (size_t n = 1; n < N - 1; n++)
    {
        x[n] = d * (bed->P[n] - bed->P[n - 1]);
    }

    x[N - 1] = bed->step->outletVelocityRHS(bed);

}

void ContinuityVelocitySystem::correctVelocity() {
    
    const double d = -bed->kappa / (bed->dx * bed->viscosity);
    const auto N = (int)x.size();

    for (int n = 1; n < N - 1; n++)
    {
        x[n] = x[n] + 0.7 * d * (bed->dP[n] - bed->dP[n - 1]);
    }

    x[N - 1] = x[N - 1] + 0.7 * 2 * d * (bed->dP[N - 2]);

}
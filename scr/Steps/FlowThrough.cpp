#include "FlowThrough.h"
#include <vector>

double FlowThrough::inletPressureRHS(const std::vector<std::vector<double>>& xPrev, const double& dt, const double& dx, const std::vector<double>& x, const std::vector<double>& u_x, const double& vis, const double& kappa) {
    return (xPrev[0][0] -
        (2 * dt / (dx)) *
        ((x[0] * u_x[0]) - (inletVelocity->boundaryValue * x[0] + dx * vis * inletVelocity->boundaryValue * inletVelocity->boundaryValue / (2 * kappa))));
}

double FlowThrough::outletPressureRHS() {
    // None required
    return 0;
}

double FlowThrough::inletVelocityRHS(const double& kappa, const double& dx, const double& vis, const std::vector<double>& p_x) {
    return (inletVelocity->boundaryValue / 2 -
        (kappa / (2 * dx * vis)) * (p_x[1] - p_x[0]));
}

double FlowThrough::outletVelocityRHS(const double& kappa, const double& dx, const double& vis, const std::vector<double>& p_x, const int& lastIndex) {
    return ((-kappa / (dx * vis)) * (outletPressure->boundaryValue -
        0.5 * (p_x[lastIndex] + p_x[lastIndex - 1])));
}

void FlowThrough::updateBoundaryConditions(const double& currentTime) {
    inletVelocity->update(currentTime);
    outletPressure->update(currentTime);
}

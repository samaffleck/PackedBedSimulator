#include "FlowThrough.h"
#include <vector>

double FlowThrough::inletPressureRHS(const std::vector<std::vector<double>>& xPrev, const double& dt, const double& dx, const std::vector<double>& x, const std::vector<double>& u_x, const std::vector<double>& t_x, const double& vis, const double& kappa) {
    return (xPrev[0][0] -
        (2 * dt / (dx)) *
        ((x[0] * u_x[0]) - (inletVelocity->boundaryValue * x[0] + dx * vis * inletVelocity->boundaryValue * inletVelocity->boundaryValue / (2 * kappa * R * t_x[0]))));
}

double FlowThrough::outletPressureRHS() {
    // None required
    return 0;
}

double FlowThrough::inletVelocityRHS(const double& kappa, const double& dx, const double& vis, const std::vector<double>& c_x, const std::vector<double>& t_x) {
    return (inletVelocity->boundaryValue / 2 -
        (kappa * R / (2 * dx * vis)) * (c_x[1] * t_x[1] - c_x[0] * t_x[0]));
}

double FlowThrough::outletVelocityRHS(const double& kappa, const double& dx, const double& vis, const std::vector<double>& c_x, const std::vector<double>& t_x, const int& lastIndex) {
    return ((-kappa / (dx * vis)) * (outletPressure->boundaryValue -
        0.5 * R * (c_x[lastIndex] * t_x[lastIndex] + c_x[lastIndex - 1] * t_x[lastIndex - 1])));
}

void FlowThrough::updateBoundaryConditions(const double& currentTime) {
    inletVelocity->update(currentTime);
    outletPressure->update(currentTime);
}

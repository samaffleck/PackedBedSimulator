#include "Pressurize.h"


double Pressurize::inletPressureRHS(const std::vector<std::vector<double>>& xPrev, const double& dt, const double& dx, const std::vector<double>& x, const std::vector<double>& u_x, const std::vector<double>& t_x, const double& vis, const double& kappa) {
    return (xPrev[0][0] -
        (dt / dx) * (0.5 * (x[0] * u_x[0] + x[1] * u_x[1]) +
            (2 * kappa * inletPressure->boundaryValue / (vis * R * t_x[0] * dx)) *
            (x[0] * t_x[0] * R - inletPressure->boundaryValue)));
}

double Pressurize::outletPressureRHS() {
    // None required
    return 0;
}

double Pressurize::inletVelocityRHS(const double& kappa, const double& dx, const double& vis, const std::vector<double>& c_x, const std::vector<double>& t_x) {
    return ((-kappa / (dx * vis)) * (0.5 * R * (c_x[0] * t_x[0] + c_x[1] * t_x[1]) -
        inletPressure->boundaryValue));
}

double Pressurize::outletVelocityRHS(const double& kappa, const double& dx, const double& vis, const std::vector<double>& c_x, const std::vector<double>& t_x, const int& lastIndex) {
    return ((-kappa * R / (2 * dx * vis)) * (c_x[lastIndex] * t_x[lastIndex] - c_x[lastIndex - 1] * t_x[lastIndex - 1]));
}

void Pressurize::updateBoundaryConditions(const double& currentTime) {
    inletPressure->update(currentTime);
}

#include "Pressurize.h"


double Pressurize::inletPressureRHS(const std::vector<std::vector<double>>& xPrev, const double& dt, const double& dx, const std::vector<double>& x, const std::vector<double>& u_x, const double& vis, const double& kappa) {
    return (xPrev[0][0] - 
            (dt / dx) * (0.5 * (x[0] * u_x[0] + x[1] * u_x[1]) -
            inletPressure->boundaryValue * ((2 * kappa) / (dx * vis)) * 
            (inletPressure->boundaryValue - x[0])));
}

double Pressurize::outletPressureRHS() {
    // None required
    return 0;
}

double Pressurize::inletVelocityRHS(const double& kappa, const double& dx, const double& vis, const std::vector<double>& p_x) {
    return ((-kappa / (dx * vis)) * (0.5 * (p_x[0] + p_x[1]) - 
        inletPressure->boundaryValue));
}

double Pressurize::outletVelocityRHS(const double& kappa, const double& dx, const double& vis, const std::vector<double>& p_x, const int& lastIndex) {
    return ((-kappa / (2 * dx * vis)) * (p_x[lastIndex] - p_x[lastIndex - 1]));
}

void Pressurize::updateBoundaryConditions(const double& currentTime) {
    inletPressure->update(currentTime);
}

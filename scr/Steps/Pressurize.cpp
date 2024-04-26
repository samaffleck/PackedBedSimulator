#include "Pressurize.h"
#include "../SystemObjects/PackedBed.h"

double Pressurize::inletDensityRHS(PackedBed* bed, const std::vector<double>& x, const std::vector<std::vector<double>>& xPrev, const double& dt) {
    return (xPrev[0][0] -
        (dt / bed->dx) * (0.5 * (x[0] * bed->U[0] + x[1] * bed->U[1]) +
            (2 * bed->kappa * inletPressure->boundaryValue / (bed->viscosity * R * bed->T[0] * bed->dx)) *
            (x[0] * bed->T[0] * R - inletPressure->boundaryValue)));
}

double Pressurize::outletDensityRHS(PackedBed* bed, const std::vector<double>& x, const std::vector<std::vector<double>>& xPrev, const double& dt) {
    // None required
    return 0;
}

double Pressurize::inletVelocityRHS(PackedBed* bed) {
    return ((-bed->kappa / (bed->dx * bed->viscosity)) * (0.5 * R * (bed->C[0] * bed->T[0] + bed->C[1] * bed->T[1]) -
        inletPressure->boundaryValue));
}

double Pressurize::outletVelocityRHS(PackedBed* bed) {
    return ((-bed->kappa * R / (2 * bed->dx * bed->viscosity)) * (bed->C[bed->numberOfCells - 1] * bed->T[bed->numberOfCells - 1] - bed->C[bed->numberOfCells - 2] * bed->T[bed->numberOfCells - 2]));
}

double Pressurize::inletTemperatureRHS(PackedBed* bed, const std::vector<double>& x, const std::vector<std::vector<double>>& xPrev, const double& dt) {
    return (xPrev[0][0] -
        (2 * dt / bed->dx) * bed->U[0] * (x[0] - inletTemperature->boundaryValue) +
        (dt / (bed->dx * bed->dx)) * bed->lambda * (x[1] - 3 * x[0] + 2 * inletTemperature->boundaryValue));
}

double Pressurize::outletTemperatureRHS(PackedBed* bed, const std::vector<double>& x, const std::vector<std::vector<double>>& xPrev, const double& dt) {
    int i = bed->numberOfCells - 1;
    return (xPrev[0][i] -
        (dt / bed->dx) * bed->U[i] * (x[i] - x[i - 1]) +
        (dt / (bed->dx * bed->dx)) * bed->lambda * (x[i - 1] - x[i]));
}

void Pressurize::updateBoundaryConditions(const double& currentTime) {
    inletPressure->update(currentTime);
}

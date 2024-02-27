#include "FlowThrough.h"
#include "../SystemObjects/PackedBed.h"
#include <vector>

double FlowThrough::inletDensityRHS(PackedBed* bed, const std::vector<double>& x, const std::vector<std::vector<double>>& xPrev, const double& dt) {
    return (xPrev[0][0] -
        (2 * dt / (bed->dx)) *
        ((x[0] * bed->U[0]) - (inletVelocity->boundaryValue * x[0] + bed->dx * bed->viscosity * inletVelocity->boundaryValue * inletVelocity->boundaryValue / (2 * bed->kappa * R * bed->T[0]))));
}

double FlowThrough::outletDensityRHS() {
    // None required
    return 0;
}

double FlowThrough::inletVelocityRHS(PackedBed* bed) {
    return (inletVelocity->boundaryValue / 2 -
        (bed->kappa * R / (2 * bed->dx * bed->viscosity)) * (bed->C[1] * bed->T[1] - bed->C[0] * bed->T[0]));
}

double FlowThrough::outletVelocityRHS(PackedBed* bed) {
    return ((-bed->kappa / (bed->dx * bed->viscosity)) * (outletPressure->boundaryValue -
        0.5 * R * (bed->C[bed->numberOfCells - 1] * bed->T[bed->numberOfCells - 1] + bed->C[bed->numberOfCells - 2] * bed->T[bed->numberOfCells - 2])));
}

double FlowThrough::inletTemperatureRHS(PackedBed* bed, const std::vector<double>& x, const std::vector<std::vector<double>>& xPrev, const double& dt) {
    return (xPrev[0][0] -
        (2 * dt / bed->dx) * bed->U[0] * (x[0] - inletTemperature->boundaryValue) +
        (dt / (bed->dx * bed->dx)) * bed->lambda * (x[1] - 3 * x[0] + 2 * inletTemperature->boundaryValue));
}

double FlowThrough::outletTemperatureRHS(PackedBed* bed, const std::vector<double>& x, const std::vector<std::vector<double>>& xPrev, const double& dt) {
    int i = bed->numberOfCells - 1;
    return (xPrev[0][i] -
        (dt / bed->dx) * bed->U[i] * (x[i] - x[i - 1]) +
        (dt / (bed->dx * bed->dx)) * bed->lambda * (x[i - 1] - x[i]));
}

void FlowThrough::updateBoundaryConditions(const double& currentTime) {
    inletVelocity->update(currentTime);
    outletPressure->update(currentTime);
}

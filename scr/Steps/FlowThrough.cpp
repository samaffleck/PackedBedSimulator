#include "FlowThrough.h"
#include "../SystemObjects/PackedBed.h"
#include <vector>

double FlowThrough::inletDensityRHS(PackedBed* bed, const std::vector<double>& x, const std::vector<std::vector<double>>& xPrev, const double& dt) {
    
    // staggered grid
    /*
    const double Tin = inletTemperature->boundaryValue;
    const double Uin = inletVelocity->boundaryValue;
    const double vis = bed->viscosity;
    const double dx = bed->dx;
    const double kap = bed->kappa;

    double Cinlet = (1 / (R * Tin)) * 
        (bed->P[0] + bed->referencePressure + (Uin * vis * dx / (2 * kap)));
    //return Cinlet * Uin / dx;
    //return bed->C[0] * Uin / dx;
    */

    // Density
    const double molarFluxIn = inletVelocity->boundaryValue; // TODO: This is molar flux, so rename to molar flux instead of inlet velocity
    return molarFluxIn / bed->dx;

}

double FlowThrough::outletDensityRHS(PackedBed* bed, const std::vector<double>& x, const std::vector<std::vector<double>>& xPrev, const double& dt) {
    
    // Stagered pressure source 
    /*
    const int N = bed->numberOfCells;
    const double dx_inv = 1 / bed->dx;
    const double ap = -bed->viscosity * bed->dx / bed->kappa;
    const double Pout = outletPressure->boundaryValue;
    const double UN = (2 / ap) * (Pout - bed->P[N - 1]);
    //const double UN = (1 / ap) * (6 * Pout - 7 * bed->P[N - 1] + bed->P[N - 2]);

    //return -(UN * dx_inv * (Pout + bed->referencePressure) / (R * bed->T[N - 1]));
    */

    // Stagered density source
    const int N = bed->numberOfCells;
    const double ap = -bed->viscosity * bed->dx / bed->kappa;
    const double Pout = outletPressure->boundaryValue;
    const double UN = bed->U[N];
    //const double UN = (2 / ap) * (Pout - bed->P[N - 1]);
    const double Cout = (Pout + bed->referencePressure) / (R * bed->T[N - 1]);

    return -(UN * Cout / bed->dx);

}

double FlowThrough::inletVelocityRHS(PackedBed* bed) {
    // For staggered grid
    //return inletVelocity->boundaryValue;

    // Density solver on stagered grid
    const double ap = bed->viscosity * bed->dx / (2 * bed->kappa);
    const double molarFluxIn = inletVelocity->boundaryValue; // TODO: rename.
    double Uin = molarFluxIn / bed->C[0]; // Initial guess;
    double Cin = (bed->P[0] + bed->referencePressure + ap * Uin ) / (R * inletTemperature->boundaryValue);
    double Uin_old = 100 * Uin; // To get in the loop;
    
    while ((Uin - Uin_old) / Uin > 1e-5) {
        Uin_old = Uin;
        Uin = molarFluxIn / Cin;
        Cin = (bed->P[0] + bed->referencePressure + ap * Uin) / (R * inletTemperature->boundaryValue);
    }

    return Uin;

}

double FlowThrough::outletVelocityRHS(PackedBed* bed) {
    // For staggered grid
    const int N = bed->numberOfCells;
    const double d = -bed->kappa / (bed->viscosity * bed->dx);
    const double Pout = outletPressure->boundaryValue;
    
    return 2 * d * (Pout - bed->P[N - 1]);
    
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

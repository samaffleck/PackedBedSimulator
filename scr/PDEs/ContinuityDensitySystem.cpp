#include "ContinuityDensitySystem.h"


void ContinuityDensitySystem::updateRHS(const double& dt) {

    const double dtInv = 1 / dt;
    const double dxInv = 1 / bed->dx;

    rhs[0] = bed->step->inletDensityRHS(bed, x, xPrev, dt);

    for (size_t i = 1; i < x.size(); i++)
    {

        rhs[i] = xPrev[0][i] - (dt / (bed->et * bed->dx)) *
            ((x[i] * bed->U[i]) - (x[i - 1] * bed->U[i - 1]));


        //rhs[i] = (xPrev[0][i] * dtInv + 
        //    x[i - 1] * bed->U[i] * 0.5 * dxInv - 
        //    x[i + 1] * bed->U[i + 1] * 0.5 * dxInv) / 
        //    (dtInv + 0.5 * dxInv * (bed->U[i + 1] - bed->U[i]));

    
        //rhs[i] = (4/3)*xPrev[0][i] - (1/3)* xPrev[1][i] - 
        //    ((2/3)*dt / (bed->et * bed->dx)) *
        //    ((x[i] * bed->U[i]) - (x[i - 1] * bed->U[i - 1]));
    
        //rhs[i] = (18 / 11) * xPrev[0][i] - (9 / 11) * xPrev[1][i] + (2 / 11) * xPrev[2][i] - ((6 / 11) * dt / (bed->et * bed->dx)) *
        //    ((x[i] * bed->U[i]) - (x[i - 1] * bed->U[i - 1]));

    }

    // stagered grid
    //rhs[x.size() - 1] = bed->step->outletDensityRHS(bed, x, xPrev, dt);

}

void ContinuityDensitySystem::updateLinkCoefficients(const double& dt)
{
    //BDF_1(dt);
    BDF_1_Density(dt);
    //BDF_2(dt);
    //BDF_3(dt);
    //BDF_4(dt);
    //BDF_5(dt);
    //BDF_6(dt);
    //CN(dt);
    //SIMPLE(dt);

}


void ContinuityDensitySystem::SIMPLE(const double& dt) {

    const auto N = (int)x.size();
    const double dx_inv = 1 / bed->dx;
    const double et_dt = bed->et / dt;
    const double alpha = et_dt / (R * bed->T[0]);
    const double ap = -bed->viscosity * bed->dx / bed->kappa;
    const double dx_ap_inv = 1 / (ap * bed->dx);
    const double relax = 1;
    double Cw = 0;
    double Ce = 0;

    // Inlet node
    Ce = 0.5 * (bed->C[0] + bed->C[1]);   // CD
    //Ce = bed->C[0];     // UDS

    // 1st Order BDF
    Aw[0] = 0;  // Not used
    //Ao[0] = alpha - Ce * dx_ap_inv;
    Ao[0] = (alpha - Ce * dx_ap_inv) / relax;
    Ae[0] = Ce * dx_ap_inv;
    So[0] = alpha * (bed->P_old[0] - bed->P[0]) - 
        Ce * bed->U[1] * dx_inv +
        bed->step->inletDensityRHS(bed, x, xPrev, dt) + 
        (1 - relax) * Ao[0] * x_LastItter[0];

    // Interior nodes
    for (int n = 1; n < N - 1; n++)
    {
        // CD
        Cw = 0.5 * (bed->C[n] + bed->C[n - 1]);
        Ce = 0.5 * (bed->C[n] + bed->C[n + 1]);

        // UD
        //Cw = bed->C[n - 1];
        //Ce = bed->C[n];

        // 1st Order BDF
        Aw[n] = Cw * dx_ap_inv;
        Ao[n] = alpha - Ce * dx_ap_inv - Cw * dx_ap_inv;
        Ae[n] = Ce * dx_ap_inv;
        So[n] = alpha * (bed->P_old[n] - bed->P[n]) + 
            Cw * bed->U[n] * dx_inv - Ce * bed->U[n + 1] * dx_inv + 
            (1 - relax) * Ao[n] * x_LastItter[n];
        
    }

    // Outlet node
    // CD
    Cw = 0.5 * (bed->C[N - 1] + bed->C[N - 2]);
    
    //UDS
    //Cw = bed->C[N - 2];

    // 1st Order BDF
    Aw[N - 1] = Cw * dx_ap_inv;
    //Ao[N - 1] = alpha - Cw * dx_ap_inv;
    Ao[N - 1] = (alpha - Cw * dx_ap_inv) / relax;
    Ae[N - 1] = 0;  // Not used
    So[N - 1] = alpha * (bed->P_old[N - 1] - bed->P[N - 1]) + 
        Cw * bed->U[N - 1] * dx_inv +
        bed->step->outletDensityRHS(bed, x, xPrev, dt) +
        (1 - relax) * Ao[N - 1] * x_LastItter[N - 1];
    
}


void ContinuityDensitySystem::BDF_1_Density(const double& dt) {

    const auto N = (int)x.size();
    const double alpha = bed->et / dt;
    const double gamma = 1.0;
    const double gamma1 = (gamma + 1) / 2;
    const double gamma2 = (1 - gamma) / 2;
    const double gamma1_dx = gamma1 / bed->dx;
    const double gamma2_dx = gamma2 / bed->dx;

    // Inlet node    
    Aw[0] = 0;  // Not used
    Ao[0] = alpha + gamma1_dx * bed->U[1];
    Ae[0] = gamma2_dx * bed->U[1];
    So[0] = xPrev[0][0] * alpha + bed->step->inletDensityRHS(bed, x, xPrev, dt);

    // Interior nodes
    for (int n = 1; n < N - 1; n++)
    {
        Aw[n] = - gamma1_dx * bed->U[n];
        Ao[n] = alpha + gamma1_dx * bed->U[n + 1] - gamma2_dx * bed->U[n];
        Ae[n] = gamma2_dx * bed->U[n + 1];
        So[n] = xPrev[0][n] * alpha;
    }

    // Outlet node
    Aw[N - 1] = - gamma1_dx * bed->U[N - 1];
    Ao[N - 1] = alpha - gamma2_dx * bed->U[N - 1];
    Ae[N - 1] = 0;  // Not used
    So[N - 1] = xPrev[0][N - 1] * alpha + bed->step->outletDensityRHS(bed, x, xPrev, dt);

}


void ContinuityDensitySystem::BDF_1(const double& dt) {

    const auto N = (int)x.size();
    const double alpha = bed->et / (R * bed->T[0] * dt);
    const double beta = bed->kappa / (bed->viscosity * bed->dx * bed->dx);
    double Cw = 0;
    double Ce = 0;

    // Inlet node
    Ce = 0.5 * (bed->C[0] + bed->C[1]);

    // 1st Order BDF
    Aw[0] = 0;  // Not used
    Ao[0] = alpha + beta * Ce;
    Ae[0] = - Ce * beta;
    So[0] = xPrev[0][0] * alpha + bed->step->inletDensityRHS(bed, x, xPrev, dt);

    // Interior nodes
    for (int n = 1; n < N - 1; n++)
    {
        Cw = 0.5 * (bed->C[n] + bed->C[n - 1]);
        Ce = 0.5 * (bed->C[n] + bed->C[n + 1]);

        // 1st Order BDF
        Aw[n] = - Cw * beta;
        Ao[n] = alpha + beta * (Cw + Ce);
        Ae[n] = - Ce * beta;
        So[n] = xPrev[0][n] * alpha;

    }

    // Outlet node
    Cw = 0.5 * (bed->C[N - 1] + bed->C[N - 2]);

    // 1st Order BDF
    Aw[N - 1] = - Cw * beta;
    Ao[N - 1] = alpha + beta * Cw;
    Ae[N - 1] = 0;  // Not used
    So[N - 1] = xPrev[0][N - 1] * alpha + bed->step->outletDensityRHS(bed, x, xPrev, dt);

}


void ContinuityDensitySystem::BDF_2(const double& dt) {

    const auto N = (int)x.size();
    const double alpha = 3 * bed->et / (2 * R * bed->T[0] * dt);
    const double beta = bed->kappa / (bed->viscosity * bed->dx * bed->dx);
    double Cw = 0;
    double Ce = 0;

    // Inlet node
    Ce = 0.5 * (bed->C[0] + bed->C[1]);

    // 1st Order BDF
    Aw[0] = 0;  // Not used
    Ao[0] = alpha + beta * Ce;
    Ae[0] = -Ce * beta;
    So[0] = alpha * ((4 / 3) * xPrev[0][0] - (1 / 3) * xPrev[1][0]) + 
        bed->step->inletDensityRHS(bed, x, xPrev, dt);

    // Interior nodes
    for (int n = 1; n < N - 1; n++)
    {
        Cw = 0.5 * (bed->C[n] + bed->C[n - 1]);
        Ce = 0.5 * (bed->C[n] + bed->C[n + 1]);

        // 1st Order BDF
        Aw[n] = -Cw * beta;
        Ao[n] = alpha + beta * (Cw + Ce);
        Ae[n] = -Ce * beta;
        So[n] = alpha * ((4 / 3) * xPrev[0][n] - (1 / 3) * xPrev[1][n]);

    }

    // Outlet node
    Cw = 0.5 * (bed->C[N - 1] + bed->C[N - 2]);

    // 1st Order BDF
    Aw[N - 1] = -Cw * beta;
    Ao[N - 1] = alpha + beta * Cw;
    Ae[N - 1] = 0;  // Not used
    So[N - 1] = alpha * ((4 / 3) * xPrev[0][N - 1] - (1 / 3) * xPrev[1][N - 1]) + 
        bed->step->outletDensityRHS(bed, x, xPrev, dt);

}


void ContinuityDensitySystem::BDF_3(const double& dt) {

    const auto N = (int)x.size();
    const double alpha = 11 * bed->et / (6 * R * bed->T[0] * dt);
    const double beta = bed->kappa / (bed->viscosity * bed->dx * bed->dx);
    double Cw = 0;
    double Ce = 0;

    // Inlet node
    Ce = 0.5 * (bed->C[0] + bed->C[1]);

    // 1st Order BDF
    Aw[0] = 0;  // Not used
    Ao[0] = alpha + beta * Ce;
    Ae[0] = -Ce * beta;
    So[0] = alpha * ((18 / 11) * xPrev[0][0] - (9 / 11) * xPrev[1][0] + (2 / 11) * xPrev[2][0]) +
        bed->step->inletDensityRHS(bed, x, xPrev, dt);

    // Interior nodes
    for (int n = 1; n < N - 1; n++)
    {
        Cw = 0.5 * (bed->C[n] + bed->C[n - 1]);
        Ce = 0.5 * (bed->C[n] + bed->C[n + 1]);

        // 1st Order BDF
        Aw[n] = -Cw * beta;
        Ao[n] = alpha + beta * (Cw + Ce);
        Ae[n] = -Ce * beta;
        So[n] = alpha * ((18 / 11) * xPrev[0][n] - (9 / 11) * xPrev[1][n] + (2 / 11) * xPrev[2][n]);

    }

    // Outlet node
    Cw = 0.5 * (bed->C[N - 1] + bed->C[N - 2]);

    // 1st Order BDF
    Aw[N - 1] = -Cw * beta;
    Ao[N - 1] = alpha + beta * Cw;
    Ae[N - 1] = 0;  // Not used
    So[N - 1] = alpha * ((18 / 11) * xPrev[0][N - 1] - (9 / 11) * xPrev[1][N - 1] + (2 / 11) * xPrev[2][N - 1]) +
        bed->step->outletDensityRHS(bed, x, xPrev, dt);

}


void ContinuityDensitySystem::BDF_4(const double& dt) {

    const auto N = (int)x.size();
    const double alpha = 25 * bed->et / (12 * R * bed->T[0] * dt);
    const double beta = bed->kappa / (bed->viscosity * bed->dx * bed->dx);
    double Cw = 0;
    double Ce = 0;

    // Inlet node
    Ce = 0.5 * (bed->C[0] + bed->C[1]);

    // 1st Order BDF
    Aw[0] = 0;  // Not used
    Ao[0] = alpha + beta * Ce;
    Ae[0] = -Ce * beta;
    So[0] = alpha * ((48 / 25) * xPrev[0][0] - (36 / 25) * xPrev[1][0] + (16 / 25) * xPrev[2][0] - (3 / 25) * xPrev[3][0]) +
        bed->step->inletDensityRHS(bed, x, xPrev, dt);

    // Interior nodes
    for (int n = 1; n < N - 1; n++)
    {
        Cw = 0.5 * (bed->C[n] + bed->C[n - 1]);
        Ce = 0.5 * (bed->C[n] + bed->C[n + 1]);

        // 1st Order BDF
        Aw[n] = -Cw * beta;
        Ao[n] = alpha + beta * (Cw + Ce);
        Ae[n] = -Ce * beta;
        So[n] = alpha * ((48 / 25) * xPrev[0][n] - (36 / 25) * xPrev[1][n] + (16 / 25) * xPrev[2][n] - (3 / 25) * xPrev[3][n]);

    }

    // Outlet node
    Cw = 0.5 * (bed->C[N - 1] + bed->C[N - 2]);

    // 1st Order BDF
    Aw[N - 1] = -Cw * beta;
    Ao[N - 1] = alpha + beta * Cw;
    Ae[N - 1] = 0;  // Not used
    So[N - 1] = alpha * ((48 / 25) * xPrev[0][N - 1] - (36 / 25) * xPrev[1][N - 1] + (16 / 25) * xPrev[2][N - 1] - (3 / 25) * xPrev[3][N - 1]) +
        bed->step->outletDensityRHS(bed, x, xPrev, dt);

}


void ContinuityDensitySystem::BDF_5(const double& dt) {

    const auto N = (int)x.size();
    const double alpha = 137 * bed->et / (60 * R * bed->T[0] * dt);
    const double beta = bed->kappa / (bed->viscosity * bed->dx * bed->dx);
    double Cw = 0;
    double Ce = 0;

    // Inlet node
    Ce = 0.5 * (bed->C[0] + bed->C[1]);

    // 1st Order BDF
    Aw[0] = 0;  // Not used
    Ao[0] = alpha + beta * Ce;
    Ae[0] = -Ce * beta;
    So[0] = alpha * ((300 / 137) * xPrev[0][0] - (300 / 137) * xPrev[1][0] + (200 / 137) * xPrev[2][0] - (75 / 137) * xPrev[3][0] + (12 / 137) * xPrev[4][0]) +
        bed->step->inletDensityRHS(bed, x, xPrev, dt);

    // Interior nodes
    for (int n = 1; n < N - 1; n++)
    {
        Cw = 0.5 * (bed->C[n] + bed->C[n - 1]);
        Ce = 0.5 * (bed->C[n] + bed->C[n + 1]);

        // 1st Order BDF
        Aw[n] = -Cw * beta;
        Ao[n] = alpha + beta * (Cw + Ce);
        Ae[n] = -Ce * beta;
        So[n] = alpha * ((300 / 137) * xPrev[0][n] - (300 / 137) * xPrev[1][n] + (200 / 137) * xPrev[2][n] - (75 / 137) * xPrev[3][n] + (12 / 137) * xPrev[4][n]);

    }

    // Outlet node
    Cw = 0.5 * (bed->C[N - 1] + bed->C[N - 2]);

    // 1st Order BDF
    Aw[N - 1] = -Cw * beta;
    Ao[N - 1] = alpha + beta * Cw;
    Ae[N - 1] = 0;  // Not used
    So[N - 1] = alpha * ((300 / 137) * xPrev[0][N - 1] - (300 / 137) * xPrev[1][N - 1] + (200 / 137) * xPrev[2][N - 1] - (75 / 137) * xPrev[3][N - 1] + (12 / 137) * xPrev[4][N - 1]) +
        bed->step->outletDensityRHS(bed, x, xPrev, dt);

}


void ContinuityDensitySystem::BDF_6(const double& dt) {

    const auto N = (int)x.size();
    const double alpha = 147 * bed->et / (60 * R * bed->T[0] * dt);
    const double beta = bed->kappa / (bed->viscosity * bed->dx * bed->dx);
    double Cw = 0;
    double Ce = 0;

    // Inlet node
    Ce = 0.5 * (bed->C[0] + bed->C[1]);

    // 1st Order BDF
    Aw[0] = 0;  // Not used
    Ao[0] = alpha + beta * Ce;
    Ae[0] = -Ce * beta;
    So[0] = alpha * ((360 / 147) * xPrev[0][0] - (450 / 147) * xPrev[1][0] + (400 / 147) * xPrev[2][0] - (225 / 147) * xPrev[3][0] + (72 / 147) * xPrev[4][0] - (10 / 147) * xPrev[5][0]) +
        bed->step->inletDensityRHS(bed, x, xPrev, dt);

    // Interior nodes
    for (int n = 1; n < N - 1; n++)
    {
        Cw = 0.5 * (bed->C[n] + bed->C[n - 1]);
        Ce = 0.5 * (bed->C[n] + bed->C[n + 1]);

        // 1st Order BDF
        Aw[n] = -Cw * beta;
        Ao[n] = alpha + beta * (Cw + Ce);
        Ae[n] = -Ce * beta;
        So[n] = alpha * ((360 / 147) * xPrev[0][n] - (450 / 147) * xPrev[1][n] + (400 / 147) * xPrev[2][n] - (225 / 147) * xPrev[3][n] + (72 / 147) * xPrev[4][n] - (10 / 147) * xPrev[5][n]);

    }

    // Outlet node
    Cw = 0.5 * (bed->C[N - 1] + bed->C[N - 2]);

    // 1st Order BDF
    Aw[N - 1] = -Cw * beta;
    Ao[N - 1] = alpha + beta * Cw;
    Ae[N - 1] = 0;  // Not used
    So[N - 1] = alpha * ((360 / 147) * xPrev[0][N - 1] - (450 / 147) * xPrev[1][N - 1] + (400 / 147) * xPrev[2][N - 1] - (225 / 147) * xPrev[3][N - 1] + (72 / 147) * xPrev[4][N - 1] - (10 / 147) * xPrev[5][N - 1]) +
        bed->step->outletDensityRHS(bed, x, xPrev, dt);

}


void ContinuityDensitySystem::CN(const double& dt) {

    const auto N = (int)x.size();
    const double half_dx_inv = 1 / (2 * bed->dx);
    const double et_dt = bed->et / dt;
    const double alpha = et_dt / (R * bed->T[0]);
    const double ap = -bed->viscosity * bed->dx / bed->kappa;
    const double dx_ap_inv = 1 / (ap * bed->dx);
    double Cw = 0;
    double Ce = 0;

    // Inlet node
    Ce = 0.5 * (bed->C[0] + bed->C[1]);

    // 2nd Order CN
    Aw[0] = 0;  // Not used
    Ao[0] = alpha - Ce * dx_ap_inv / 2;
    Ae[0] = Ce * dx_ap_inv / 2;
    So[0] = xPrev[0][0] * alpha + 
        Ce * dx_ap_inv * xPrev[0][0] / 2 -
        Ce * dx_ap_inv * xPrev[0][1] / 2 +
        bed->step->inletDensityRHS(bed, x, xPrev, dt);


    // Interior nodes
    for (int n = 1; n < N - 1; n++)
    {
        Cw = 0.5 * (bed->C[n] + bed->C[n - 1]);
        Ce = 0.5 * (bed->C[n] + bed->C[n + 1]);

        // 2nd Order CN
        Aw[n] = Cw * dx_ap_inv / 2;
        Ao[n] = alpha - Ce * dx_ap_inv / 2 - Cw * dx_ap_inv / 2;
        Ae[n] = Ce * dx_ap_inv / 2;
        So[n] = xPrev[0][n] * alpha - 
            Ce * dx_ap_inv * xPrev[0][n + 1] / 2 +
            Ce * dx_ap_inv * xPrev[0][n] / 2 + 
            Cw * dx_ap_inv * xPrev[0][n] / 2 - 
            Cw * dx_ap_inv * xPrev[0][n - 1] / 2;

    }

    // Outlet node
    Cw = 0.5 * (bed->C[N - 1] + bed->C[N - 2]);

    // 2nd Order CN
    Aw[N - 1] = Cw * dx_ap_inv / 2;
    Ao[N - 1] = alpha - Cw * dx_ap_inv / 2;
    Ae[N - 1] = 0;  // Not used
    So[N - 1] = xPrev[0][N - 1] * alpha -
        Cw * dx_ap_inv * xPrev[0][N - 2] + 
        Cw * dx_ap_inv * xPrev[0][N - 1] +
        bed->step->outletDensityRHS(bed, x, xPrev, dt);

}


void ContinuityDensitySystem::correctPressure() {

    for (int n = 0; n < bed->numberOfCells; n++) {
        bed->P[n] = bed->P[n] + 0.7 * x[n];
    }

}

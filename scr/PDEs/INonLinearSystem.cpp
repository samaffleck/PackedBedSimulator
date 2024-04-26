#include "INonLinearSystem.h"
#include "../Solver/ThomasAlgoritms.h"

void INonLinearSystem::gaussSeidel(const double& alpha = 1) {

    const auto N = (int)x.size();
    const double damp = 0;

    // Inlet node
    x[0] = (1 - alpha) * x[0] + alpha * ((So[0] - x[1] * Ae[0]) / (Ao[0] + damp)) ;

    // Interior nodes
    for (size_t n = 1; n < N - 1; n++)
    {
        x[n] = (1 - alpha) * x[n] + alpha * ((So[n] - x[n + 1] * Ae[n] - x[n - 1] * Aw[n]) / (Ao[n] + damp));
    }

    // Outlet node
    x[N - 1] = (1 - alpha) * x[N - 1] + alpha * ((So[N - 1] - x[N - 2] * Aw[N - 1]) / (Ao[N - 1] + damp));

}


double INonLinearSystem::evaluateError() {

    const auto N = (int)x.size();

    /*
    // Inlet node
    e[0] = x[0] * Ao[0] - (So[0] - x[1] * Ae[0]);

    // Interior nodes
    for (size_t n = 1; n < N - 1; n++)
    {
        e[n] = x[n] * Ao[n] - (So[n] - x[n + 1] * Ae[n] - x[n - 1] * Aw[n]);
    }

    // Outlet node
    e[N - 1] = x[N - 1] * Ao[N - 1] - (So[N - 1] - x[N - 2] * Aw[N - 1]);
    */

    for (size_t i = 0; i < e.size(); i++)
    {        
        e[i] = (x[i] - x_LastItter[i]);
    }
    
    const double error = VectorNorm::L2Norm(e) / VectorNorm::L2Norm(x);
    //if (!std::isfinite(error)) dampingFactor /= 10;
    return error;

}


void INonLinearSystem::innerItteration(const int& maxItterations, const double& tolerance, const double& timeStep) {
    error = tolerance * 100;
    itterations = 0;

    x_LastItter = x;

    TDMA.solve(Ae, Ao, Aw, So, x);
    
    // Dampen the solution
    for (size_t i = 0; i < x.size(); i++)
    {
        x[i] = x_LastItter[i] + (dampingFactor / timeStep) * (x[i] - x_LastItter[i]);
    }


    // GS
    /*
    while (error > tolerance && itterations < maxItterations) {
        gaussSeidel(1.2);
        error = evaluateError() / bed->numberOfCells;
        x_LastItter = x;
        itterations++;
    }  
    */

}

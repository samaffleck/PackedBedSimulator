#include <iostream>
#include <vector>
#include <math.h>

#include "PDEs/ContinuityPressureSystem.h"
#include "PDEs/ContinuityVelocitySystem.h"
#include "PDEs/DiffusionSystem.h"

#include "Solver/NonLinearSolver.h"


int main(){

    int systemSize = 20;
    int maxItterations = 100;
    int order = 1;
    double tolerance = 1e-10;

    // Initial conditions
    std::vector<double> p0(systemSize);
    std::vector<double> u0(systemSize);
    std::vector<double> T0(systemSize);

    for (size_t i = 0; i < p0.size(); i++)
    {
        p0[i] = 1e5;
        u0[i] = 0.1;
        T0[i] = 200;
    }
    
    DiffusionSystem diffusionSystem(T0, order);
    
    NonLinearSolver solver(maxItterations, tolerance, tolerance, 1000, 1);

    solver.addNonLinearSystem(&diffusionSystem);
    //solver.addNonLinearSystem(&velocitySystem);
    //solver.addNonLinearSystem(&pressureSystem);

    solver.run();

    return 0;
}


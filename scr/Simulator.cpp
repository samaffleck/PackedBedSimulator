#include <vector>

#include "PDEs/ContinuityPressureSystem.h"
#include "PDEs/ContinuityVelocitySystem.h"
#include "PDEs/DiffusionSystem.h"

#include "Solver/NonLinearSolver.h"

#include "BoundaryConditions/BoundaryCondition_Constant.h"
#include "BoundaryConditions/BoundaryCondition_Flux.h"


int main(){

    int systemSize = 10;
    int maxItterations = 1000;
    int order = 1;
    double tolerance = 1e-6;

    // Initial conditions
    std::vector<double> p0(systemSize);
    std::vector<double> u0(systemSize);
    std::vector<double> T0(systemSize);

    for (size_t i = 0; i < p0.size(); i++)
    {
        p0[i] = 1e5;
        u0[i] = 0.01;
        T0[i] = 100;
    }
    
    DiffusionSystem diffusionSystem(T0, order);
    BoundaryCondition_Constant inletBoundaryCondition(100);
    BoundaryCondition_Constant outletBoundaryCondition(200);
    diffusionSystem.inletBoundaryCondition = &inletBoundaryCondition;
    diffusionSystem.outletBoundaryCondition = &outletBoundaryCondition;

    //ContinuityPressureSystem pressureSystem(p0, order);
    //BoundaryCondition_Constant outletPressure(1e5);
    //BoundaryCondition_Constant inletPressure(1.1e5);
    //pressureSystem.inletBoundaryCondition = &inletPressure;
    //pressureSystem.outletBoundaryCondition = &outletPressure;

    //ContinuityVelocitySystem velocitySystem(u0, order);
    //BoundaryCondition_Constant inletVelocity(0.5);
    //velocitySystem.inletBoundaryCondition = &inletVelocity;

    // Coupple the pressure and velocity systems
    //pressureSystem.uSystem = &velocitySystem;
    //velocitySystem.pSystem = &pressureSystem;
    
    // Set up the solver
    NonLinearSolver solver(maxItterations, tolerance, tolerance, 100, 0.1);

    solver.addNonLinearSystem(&diffusionSystem);
    //solver.addNonLinearSystem(&pressureSystem);
    //solver.addNonLinearSystem(&velocitySystem);
    
    solver.run();

    return 0;
}

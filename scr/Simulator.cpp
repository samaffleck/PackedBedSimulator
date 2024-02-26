#include <vector>

#include "PDEs/ContinuityPressureSystem.h"
#include "PDEs/ContinuityVelocitySystem.h"
#include "PDEs/DiffusionSystem.h"
#include "PDEs/BurgesEquation.h"
#include "PDEs/ConvectionEquation.h"

#include "Solver/NonLinearSolver.h"

#include "BoundaryConditions/BoundaryCondition_Constant.h"
#include "BoundaryConditions/BoundaryCondition_Flux.h"
#include "BoundaryConditions/BoundaryCondition_Step.h"
#include "BoundaryConditions/BoundaryCondition_Pulse.h"


int main(){

    int systemSize = 100;
    double length = 35; //[m]

    int maxItterations = 100;
    int order = 1;
    double tolerance = 1e-6;

    // Initial conditions
    std::vector<double> p0(systemSize);
    std::vector<double> u0(systemSize);
    std::vector<double> T0(systemSize);

    for (size_t i = 0; i < p0.size(); i++)
    {
        p0[i] = 1e5;
        u0[i] = 0.0;
        T0[i] = 100;
    }
    
    DiffusionSystem diffusionSystem(T0, order, systemSize, length);
    BoundaryCondition_Constant inletBoundaryCondition(100);
    BoundaryCondition_Constant outletBoundaryCondition(200);
    diffusionSystem.inletBoundaryCondition = &inletBoundaryCondition;
    diffusionSystem.outletBoundaryCondition = &outletBoundaryCondition;

    ConvectionEquation convectionEquation(u0, order, systemSize, length);
    BoundaryCondition_Pulse inletPulse(0,1,0,0.5);
    //BoundaryCondition_Step inletStep(0,1,0);
    convectionEquation.inletBoundaryCondition = &inletPulse;

    //BurgesEquation burgesEquation(u0, order, systemSize, length);
    //BoundaryCondition_Constant inletVelocity(0);
    //burgesEquation.inletBoundaryCondition = &inletVelocity;

    //ContinuityPressureSystem pressureSystem(p0, order, systemSize, length);
    //BoundaryCondition_Constant outletPressure(1e5);
    //BoundaryCondition_Constant inletPressure(1.1e5);
    //pressureSystem.inletBoundaryCondition = &inletPressure;
    //pressureSystem.outletBoundaryCondition = &outletPressure;

    //ContinuityVelocitySystem velocitySystem(u0, order, systemSize, length);
    //BoundaryCondition_Constant inletVelocity(0.5);
    //velocitySystem.inletBoundaryCondition = &inletVelocity;

    // Coupple the pressure and velocity systems
    //pressureSystem.uSystem = &velocitySystem;
    //velocitySystem.pSystem = &pressureSystem;
    
    // Set up the solver
    NonLinearSolver solver(maxItterations, tolerance, tolerance, 5, 0.001);

    convectionEquation.time = &solver.currentTime;

    //solver.addNonLinearSystem(&diffusionSystem);
    solver.addNonLinearSystem(&convectionEquation);
    //solver.addNonLinearSystem(&burgesEquation);
    //solver.addNonLinearSystem(&pressureSystem);
    //solver.addNonLinearSystem(&velocitySystem);
    
    solver.run();

    return 0;
}

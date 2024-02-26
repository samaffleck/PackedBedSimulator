#include <vector>

#include "SystemObjects/PackedBed.h"

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

    int maxItterations = 50;
    double tolerance = 1e-6;

    PackedBed bed(35, 0.4, 0.3, 0.001, 2e-5);
    bed.initialise(100, 298, 0, 1e5);

    BoundaryCondition_Constant outletPressure(1e5);
    BoundaryCondition_Constant inletPressure(1.1e5);
    bed.pressureSystem.inletBoundaryCondition = &inletPressure;
    bed.pressureSystem.outletBoundaryCondition = &outletPressure;

    BoundaryCondition_Constant inletVelocity(0.5);
    bed.velocitySystem.inletBoundaryCondition = &inletVelocity;


    // Initial conditions
    //std::vector<double> p0(systemSize);
    //std::vector<double> u0(systemSize);
    //std::vector<double> T0(systemSize);

    //for (size_t i = 0; i < p0.size(); i++)
    //{
    //    p0[i] = 1e5;
    //    u0[i] = 0.5;
    //    T0[i] = 100;
    //}
    
    //DiffusionSystem diffusionSystem(T0, order, systemSize, length);
    //BoundaryCondition_Constant inletBoundaryCondition(100);
    //BoundaryCondition_Constant outletBoundaryCondition(200);
    //diffusionSystem.inletBoundaryCondition = &inletBoundaryCondition;
    //diffusionSystem.outletBoundaryCondition = &outletBoundaryCondition;

    //ConvectionEquation convectionEquation(u0, order, systemSize, length);
    //BoundaryCondition_Pulse inletPulse(0,1,0,0.5);
    //BoundaryCondition_Step inletStep(0,1,0);
    //convectionEquation.inletBoundaryCondition = &inletPulse;

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
    NonLinearSolver solver(bed, maxItterations, tolerance, tolerance, 200, 0.001);

    //solver.addNonLinearSystem(&diffusionSystem);
    //solver.addNonLinearSystem(&convectionEquation);
    //solver.addNonLinearSystem(&burgesEquation);
    
    //solver.addNonLinearSystem(&bed.pressureSystem);
    //solver.addNonLinearSystem(&bed.velocitySystem);
    
    solver.run();

    return 0;
}

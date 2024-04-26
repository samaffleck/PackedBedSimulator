#include <vector>

#include "SystemObjects/PackedBed.h"
#include "Solver/NonLinearSolver.h"

#include "Steps/FlowThrough.h"
#include "Steps/Pressurize.h"

#include "BoundaryConditions/BoundaryCondition_Constant.h"
#include "BoundaryConditions/BoundaryCondition_Flux.h"
#include "BoundaryConditions/BoundaryCondition_Step.h"
#include "BoundaryConditions/BoundaryCondition_Pulse.h"
#include "BoundaryConditions/BoundaryCondition_Ramp.h"


int main(){

    int maxItterations = 500;
    double innerTolerance = 1e-2;
    double outerTolerance = 1e-6;

    PackedBed bed(0.035, 0.4, 0.3, 0.002, 2e-5);
    bed.initialise(40, 298, 0.0, 0);
    
    BoundaryCondition_Constant outletPressure(0);

    BoundaryCondition_Constant inletPressure(1e1);

    BoundaryCondition_Constant inletTemperature(298);
    
    BoundaryCondition_Constant inletVelocity(20);

    BoundaryCondition_Ramp inletPressureRamp(1e5, 1.5e5, 50);

    FlowThrough flowStep(&inletVelocity, &outletPressure, &inletTemperature);
    Pressurize pressurizeStep(&inletPressureRamp, &inletTemperature);
    
    //bed.selectStep(&pressurizeStep);
    bed.selectStep(&flowStep);

    // Set up the solver
    NonLinearSolver solver(bed, maxItterations, innerTolerance, outerTolerance, 600, 0.01); // 1e-3 to 5e-4
    solver.run();

    return 0;

}

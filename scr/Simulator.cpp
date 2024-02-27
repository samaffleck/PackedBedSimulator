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

    int maxItterations = 50;
    double innerTolerance = 1e-3;
    double outerTolerance = 1e-6;

    PackedBed bed(35, 0.4, 0.3, 0.001, 2e-5);
    bed.initialise(20, 298, 0, 1e5);
    
    BoundaryCondition_Constant outletPressure(1e5);
    BoundaryCondition_Constant inletPressure(1.1e5);
    BoundaryCondition_Constant inletVelocity(0.1);

    BoundaryCondition_Ramp inletPressureRamp(1e5, 1.5e5, 50);

    FlowThrough flowStep(&inletVelocity, &outletPressure);
    Pressurize pressurizeStep(&inletPressureRamp);
    bed.selectStep(&pressurizeStep);
    //bed.selectStep(&flowStep);

    // Set up the solver
    NonLinearSolver solver(bed, maxItterations, innerTolerance, outerTolerance, 100, 0.001);
    solver.run();

    return 0;
}

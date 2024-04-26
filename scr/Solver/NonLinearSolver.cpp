#include "NonLinearSolver.h"


void NonLinearSolver::run() {
    // Initialisation
    currentTime = 0;
    timeStep = initialTimeStep;
    dataLogger.open(3, bed.numberOfCells, bed.dx);
    
    for (int i = 0; i <= totalTime; ++i)
    {
        timeQ.push(i);
    }

    // Initialise system of equations
    //densitySystem = new ContinuityDensitySystem(&bed, bed.dP, bed.order);     // Pressure correction algo
    //densitySystem = new ContinuityDensitySystem(&bed, bed.P, bed.order);      // Direct pressure algo
    densitySystem = new ContinuityDensitySystem(&bed, bed.C, bed.order);        // Direct density algo
    velocitySystem = new ContinuityVelocitySystem(&bed, bed.U, bed.order);      
    temperatureSystem = new AdvectionDiffusionSystem(&bed, bed.T, bed.order);

    while (currentTime < totalTime)
    {
        successfulStep = false;

        // Update boundary values
        bed.step->updateBoundaryConditions(currentTime);
        
        while (successfulStep == false && timeStep > minimumTimeStep)
        {
            successfulStep = outerItteration();
        }

        if (timeStep < minimumTimeStep) break;

    }

    dataLogger.close();

    delete densitySystem;
    delete velocitySystem;
    delete temperatureSystem;

}


bool NonLinearSolver::outerItteration() {

    outerItterations = 0;
    outerError = outerTolerance * 100;

    bool isConverged = false;

    while (isConverged == false && outerItterations < maxItterations) {

        //SIMPLE();
        
        //LAPLACE();
        
        DENSITY();

        isConverged = updateSystemError();
        outerItterations++;
    }
    
    if (outerItterations == maxItterations || std::isfinite(outerError) == false) {
        rejectStep();
        return false;
    }
    else {
        acceptStep();
        return true;
    }
    
}


void NonLinearSolver::SIMPLE() {

    // solve x-momentum 
    velocitySystem->updateVelocity();

    // Set up pressure correction coeff.
    densitySystem->updateLinkCoefficients(timeStep);

    // Solve pressure correction
    densitySystem->innerItteration(maxItterations, innerTolerance, timeStep);

    // Correct pressure
    densitySystem->correctPressure();

    // Correct velocity
    velocitySystem->correctVelocity();

    // Update density
    bed.updateConstants();

}

void NonLinearSolver::LAPLACE() {

    // Set up pressure correction coeff.
    densitySystem->updateLinkCoefficients(timeStep);

    // Solve pressure correction
    densitySystem->innerItteration(maxItterations, innerTolerance, timeStep);

    // Update density
    bed.updateConstants();

}

void NonLinearSolver::DENSITY() {

    // Set up pressure correction coeff.
    densitySystem->updateLinkCoefficients(timeStep);

    // Solve pressure correction
    densitySystem->innerItteration(maxItterations, innerTolerance, timeStep);    
    
    // Update density
    bed.updateConstants();

    // Update velocity
    velocitySystem->updateVelocity();

}

void NonLinearSolver::acceptStep() {

    log(currentTime);
    currentTime += timeStep;
    //timeStep *= 2;
    if (timeStep > maxTimeStep) timeStep = maxTimeStep;
    
    // TODO: move this
    // Update the old vectors depending on the order
    for (int j = bed.order - 1; j > 0; j--)
    {
        densitySystem->xPrev[j] = densitySystem->xPrev[j - 1];
        velocitySystem->xPrev[j] = velocitySystem->xPrev[j - 1];
        temperatureSystem->xPrev[j] = temperatureSystem->xPrev[j - 1];
    }
    densitySystem->xPrev[0] = densitySystem->x;
    velocitySystem->xPrev[0] = velocitySystem->x;
    temperatureSystem->xPrev[0] = temperatureSystem->x;

    bed.P_old = bed.P;

}


void NonLinearSolver::rejectStep() {
    densitySystem->x = densitySystem->xPrev[0];
    densitySystem->x_LastItter = densitySystem->xPrev[0];
    
    velocitySystem->x = velocitySystem->xPrev[0];
    velocitySystem->x_LastItter = velocitySystem->xPrev[0];
    
    temperatureSystem->x = temperatureSystem->xPrev[0];
    temperatureSystem->x_LastItter = temperatureSystem->xPrev[0];

    bed.P = bed.P_old;

    //timeStep *= 0.5;
}


bool NonLinearSolver::updateSystemError() {

    // With all individual PDEs solved, we have to update the x vector for each system in order to get the couppled error
    //densitySystem->updateLinkCoefficients(timeStep);
    //densitySystem->gaussSeidel(1);
    //velocitySystem->gaussSeidel(timeStep);
    //temperatureSystem->gaussSeidel(timeStep);

    // Update the overall error. This must only be called after all function inner itterations are complete
    double densityNorm = densitySystem->evaluateError() / bed.numberOfCells;
    double velocityNorm = velocitySystem->evaluateError() / bed.numberOfCells;
    double temperatureNorm = temperatureSystem->evaluateError() / bed.numberOfCells;

    outerError = densityNorm;

    //if (densityNorm <= outerTolerance && velocityNorm <= outerTolerance && temperatureNorm <= outerTolerance) return true;
    //if (densityNorm <= outerTolerance && velocityNorm <= outerTolerance) return true;
    if (densityNorm <= outerTolerance) return true;
    
    return false;
}


void NonLinearSolver::log(const double& time) {
    // Log successful step
    if (time >= timeQ.front()) {
        timeQ.pop();
        std::cout << "time: " << currentTime << "\ttime step: " << timeStep << "\n";
        std::cout << "itterations: " << outerItterations << "\n";
        dataLogger.log(0, currentTime, bed.P);
        dataLogger.log(1, currentTime, bed.U);
        dataLogger.log(2, currentTime, bed.C);
    }
}

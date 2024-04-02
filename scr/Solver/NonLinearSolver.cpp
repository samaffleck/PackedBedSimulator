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
    densitySystem = new ContinuityDensitySystem(&bed, bed.C, 1);
    velocitySystem = new ContinuityVelocitySystem(&bed, bed.U, 1);
    temperatureSystem = new AdvectionDiffusionSystem(&bed, bed.T, 1);

    while (currentTime < totalTime)
    {
        successfulStep = false;

        // Update boundary values
        bed.step->updateBoundaryConditions(currentTime);
        
        while (successfulStep == false and timeStep > minimumTimeStep)
        {
            //successfulStep = outerItteration.solve(timeStep, currentTime);
            successfulStep = outerItteration();
        }

    }

    dataLogger.close();

}


bool NonLinearSolver::outerItteration() {

    outerItterations = 0;
    outerError = outerTolerance * 100;

    while (outerError > outerTolerance && outerItterations < maxItterations) {
        innerItterations();
        updateSystemError();
        outerItterations++;
    }

    if (outerItterations == maxItterations or std::isfinite(outerError) == false) {
        rejectStep();
        return false;
    }
    else {
        acceptStep();
        return true;
    }

}


void NonLinearSolver::acceptStep() {

    log(currentTime);
    currentTime += timeStep;
    timeStep *= 2;
    if (timeStep > 1) // Max time step
    {
        timeStep = 1;
    }

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

}


void NonLinearSolver::rejectStep() {
    densitySystem->x = densitySystem->xPrev[0];
    velocitySystem->x = velocitySystem->xPrev[0];
    temperatureSystem->x = temperatureSystem->xPrev[0];
    timeStep *= 0.5;
}


void NonLinearSolver::updateSystemError() {

    // Update the overall error. This must only be called after all function inner itterations are complete
    outerError = 0;
    outerError += densitySystem->evaluateError(timeStep);
    outerError += velocitySystem->evaluateError(timeStep);
    outerError += temperatureSystem->evaluateError(timeStep);
    outerError = outerError / (3 * bed.numberOfCells); // 3 is the number of equations (this may change)

}


void NonLinearSolver::innerItterations() {

    densitySystem->innerItteration(maxItterations, innerTolerance, timeStep);
    velocitySystem->innerItteration(maxItterations, innerTolerance, timeStep);
    temperatureSystem->innerItteration(maxItterations, innerTolerance, timeStep);

}


void NonLinearSolver::log(const double& time) {
    // Log successful step
    if (time >= timeQ.front()) {
        timeQ.pop();
        std::cout << "time: " << currentTime << "\ttime step: " << timeStep << "\n";
        dataLogger.log(0, currentTime, densitySystem->x);
        dataLogger.log(1, currentTime, velocitySystem->x);
        dataLogger.log(2, currentTime, temperatureSystem->x);
    }
}
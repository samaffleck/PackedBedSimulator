#pragma once

#include "../PDEs/INonLinearSystem.h"
#include "../VectorFunctions/VectorPrinter.h"


class OuterItteration {
public:
    OuterItteration(
        PackedBed& _bed,
        int _maxItterations,
        double _innerTolerance,
        double _outerTolerance) :
        bed(_bed),
        maxItterations(_maxItterations),
        innerTolerance(_innerTolerance),
        outerTolerance(_outerTolerance) {}
    ~OuterItteration() {}

    double error{};
    double innerTolerance{};
    double outerTolerance{};
    int itterations{};
    int maxItterations{};
    
    PackedBed& bed;   

    bool solve(double& timeStep, double& currentTime) {

        itterations = 0;
        error = outerTolerance * 100;

        while (error > outerTolerance && itterations < maxItterations) {
            innerItterations(timeStep);
            updateSystemError(timeStep);
            itterations++;
        }

        if (itterations == maxItterations or std::isfinite(error) == false) {
            rejectStep(timeStep);
            return false;
        }
        else {
            acceptStep(timeStep, currentTime);
            return true;
        }

    }

    void acceptStep(double& timeStep, double& currentTime) {

        currentTime += timeStep;
        timeStep *= 2;
        if (timeStep > 1) // Max time step
        {
            timeStep = 1;
        }

        // Update the old vectors depending on the order
        for (int j = bed.order - 1; j > 0; j--)
        {
            bed.pressureSystem.xPrev[j] = bed.pressureSystem.xPrev[j - 1];
            bed.velocitySystem.xPrev[j] = bed.velocitySystem.xPrev[j - 1];
        }
        bed.pressureSystem.xPrev[0] = bed.pressureSystem.x;
        bed.velocitySystem.xPrev[0] = bed.velocitySystem.x;
        
    }

    void rejectStep(double& timeStep) {
        bed.pressureSystem.x = bed.pressureSystem.xPrev[0];
        bed.velocitySystem.x = bed.velocitySystem.xPrev[0];
        timeStep *= 0.5;
    }

    void updateSystemError(double& timeStep) {
        // Update the overall error. This must only be called after all function inner itterations are complete
        error = 0;
        
        error += bed.pressureSystem.evaluateError(timeStep);
        error += bed.velocitySystem.evaluateError(timeStep);

        error = error / (2 * bed.numberOfCells); // 2 is the number of equations (this may change)

    }

    void innerItterations(double& timeStep) {
        bed.pressureSystem.innerItteration(maxItterations, innerTolerance, timeStep);
        bed.velocitySystem.innerItteration(maxItterations, innerTolerance, timeStep);
    }

};

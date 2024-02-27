#pragma once

#include "OuterItteration.h"
#include "../PDEs/INonLinearSystem.h"
#include "../DataLogger/DataLogger.h"
#include "../SystemObjects/PackedBed.h"


class NonLinearSolver {
public:
    NonLinearSolver(
        PackedBed& _bed,
        int _maxItterations,
        double _innerTolerance,
        double _outerTolerance,
        double _totalTime,
        double _initialTimeStep) :
        bed(_bed),
        outerItteration(_bed, _maxItterations, _innerTolerance, _outerTolerance),
        totalTime(_totalTime),
        initialTimeStep(_initialTimeStep) {}
    ~NonLinearSolver() {}

    double currentTime{};
    double totalTime;
    double initialTimeStep;
    double timeStep{};
    double minimumTimeStep = 0.00000001;
    bool successfulStep = false;

    PackedBed& bed;
    OuterItteration outerItteration;
    DataLogger dataLogger;

    void run() {
        // Initialisation
        currentTime = 0;
        timeStep = initialTimeStep;
        dataLogger.open(2, bed.numberOfCells, bed.dx);

        while (currentTime < totalTime)
        {
            successfulStep = false;

            // Update boundary values
            bed.densitySystem.step->updateBoundaryConditions(currentTime);
            
            while (successfulStep == false and timeStep > minimumTimeStep)
            {
                successfulStep = outerItteration.solve(timeStep, currentTime);
            }

            if (successfulStep)
            {
                // Log successful step
                std::cout << "time: " << currentTime << "\ttime step: " << timeStep << "\n";

                dataLogger.log(0, currentTime, bed.densitySystem.x);
                dataLogger.log(1, currentTime, bed.velocitySystem.x);

            }
            
        }

        dataLogger.close();

    }

};

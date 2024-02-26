#pragma once

#include "OuterItteration.h"
#include "../PDEs/INonLinearSystem.h"
#include "../DataLogger/DataLogger.h"


class NonLinearSolver {
public:
    NonLinearSolver(
        int _maxItterations,
        double _innerTolerance,
        double _outerTolerance,
        double _totalTime,
        double _initialTimeStep) :
        outerItteration(_maxItterations, _innerTolerance, _outerTolerance),
        totalTime(_totalTime),
        initialTimeStep(_initialTimeStep) {}
    ~NonLinearSolver() {}

    double currentTime{};
    double totalTime;
    double initialTimeStep;
    double timeStep{};
    double minimumTimeStep = 0.000001;
    bool successfulStep = false;

    OuterItteration outerItteration;
    DataLogger dataLogger;

    void run() {
        // Initialisation
        currentTime = 0;
        timeStep = initialTimeStep;
        dataLogger.open(outerItteration.vectorOfNonLinearSystems.size(), outerItteration.vectorOfNonLinearSystems[0]->x.size(), outerItteration.vectorOfNonLinearSystems[0]->dx);

        while (currentTime < totalTime)
        {
            successfulStep = false;

            // Update boundary values
            for (size_t i = 0; i < outerItteration.vectorOfNonLinearSystems.size(); i++)
            {
                outerItteration.vectorOfNonLinearSystems[i]->inletBoundaryCondition->update(currentTime);
            }

            while (successfulStep == false and timeStep > minimumTimeStep)
            {
                successfulStep = outerItteration.solve(timeStep, currentTime);
            }

            if (successfulStep)
            {
                currentTime += timeStep;
                //timeStep *= 2;
                if (timeStep > 1) // Max time step
                {
                    timeStep = 1;
                }

                // Update the old vectors depending on the order
                for (size_t i = 0; i < outerItteration.vectorOfNonLinearSystems.size(); i++)
                {
                    for (int j = outerItteration.vectorOfNonLinearSystems[i]->order - 1; j > 0; j--)
                    {
                        outerItteration.vectorOfNonLinearSystems[i]->xPrev[j] = outerItteration.vectorOfNonLinearSystems[i]->xPrev[j - 1];
                    }
                    outerItteration.vectorOfNonLinearSystems[i]->xPrev[0] = outerItteration.vectorOfNonLinearSystems[i]->x;
                }

                std::cout << "time: " << currentTime << "\ttime step: " << timeStep << "\n";

                for (size_t i = 0; i < outerItteration.vectorOfNonLinearSystems.size(); i++)
                {
                    dataLogger.log(i, currentTime, outerItteration.vectorOfNonLinearSystems[i]->x);
                }

            }
            else {
                /* Halve the time step and don't update x */
                //timeStep *= 0.5;
            }
            
        }

        dataLogger.close();

    }

    void addNonLinearSystem(INonLinearSystem* nonLinearSystem) {
        outerItteration.vectorOfNonLinearSystems.push_back(nonLinearSystem);
    }

};

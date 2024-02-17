#pragma once

#include "OuterItteration.h"
#include "../PDEs/INonLinearSystem.h"


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

    double currentTime;
    double totalTime;
    double initialTimeStep;
    double timeStep{};

    OuterItteration outerItteration;

    void run() {
        // Initialisation
        currentTime = 0;
        timeStep = initialTimeStep;

        while (currentTime < totalTime)
        {
            outerItteration.solve(timeStep, currentTime);

            std::cout << "time: " << currentTime << "\ttime step: " << timeStep << "\n";
        }

    }

    void addNonLinearSystem(INonLinearSystem* nonLinearSystem) {
        outerItteration.vectorOfNonLinearSystems.push_back(nonLinearSystem);
    }

};

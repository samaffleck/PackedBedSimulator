#pragma once

#include "../PDEs/INonLinearSystem.h"
#include "../DataLogger/DataLogger.h"
#include "../SystemObjects/PackedBed.h"

#include "../PDEs/AdvectionDiffusionSystem.h"
#include "../PDEs/ContinuityDensitySystem.h"
#include "../PDEs/ContinuityVelocitySystem.h"


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
        outerTolerance(_outerTolerance),
        innerTolerance(_innerTolerance),
        maxItterations(_maxItterations),
        totalTime(_totalTime),
        initialTimeStep(_initialTimeStep) {}
    ~NonLinearSolver() {
        //delete densitySystem;
        //delete velocitySystem;
        //delete temperatureSystem;
    }

    double currentTime{};
    double totalTime;
    double initialTimeStep;
    double timeStep{};
    double minimumTimeStep = 0.00000001;
    bool successfulStep = false;

    int outerItterations{};
    int maxItterations{};
    double outerError{};
    double outerTolerance{};
    double innerTolerance{};

    PackedBed& bed;
    DataLogger dataLogger;

    ContinuityDensitySystem* densitySystem = nullptr;
    ContinuityVelocitySystem* velocitySystem = nullptr;
    AdvectionDiffusionSystem* temperatureSystem = nullptr;

    void run();
    bool outerItteration();
    void acceptStep();
    void rejectStep();
    void updateSystemError();
    void innerItterations();

};

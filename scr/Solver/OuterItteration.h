#pragma once

#include "../PDEs/INonLinearSystem.h"
#include "../VectorFunctions/VectorPrinter.h"


class OuterItteration {
public:
    OuterItteration(
        int _maxItterations,
        double _innerTolerance,
        double _outerTolerance) :
        maxItterations(_maxItterations),
        innerTolerance(_innerTolerance),
        outerTolerance(_outerTolerance) {}
    ~OuterItteration() {}

    double error{};
    double innerTolerance{};
    double outerTolerance{};
    int itterations{};
    int maxItterations{};
    std::vector<INonLinearSystem*> vectorOfNonLinearSystems{};

    bool solve(double& timeStep, double& currentTime) {

        itterations = 0;
        error = outerTolerance * 100;

        while (error > outerTolerance && itterations < maxItterations) {

            for (size_t i = 0; i < vectorOfNonLinearSystems.size(); i++)
            {
                vectorOfNonLinearSystems[i]->innerItteration(maxItterations, innerTolerance, timeStep);
            }

            // Update the overall error. This must only be called after all function inner itterations are complete
            error = 0;
            for (size_t i = 0; i < vectorOfNonLinearSystems.size(); i++)
            {
                error += vectorOfNonLinearSystems[i]->evaluateError(timeStep);
            }
            error = error / (vectorOfNonLinearSystems.size() * vectorOfNonLinearSystems[0]->x.size());

            itterations++;

        }

        if (itterations == maxItterations)
        {
            /* Halve the time step and don't update x */

            //timeStep *= 0.5;

            return false; // Not a successful step

        }
        else {
            /* Accept the step, increase the time step and update x*/

            return true; // Accept the step
            //std::cout << "error = " << error << "\tOutter Itterations = " << itterations << "\n";

        }

    }

};

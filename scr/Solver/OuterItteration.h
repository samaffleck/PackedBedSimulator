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

    void solve(double& timeStep, double& currentTime) {

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

            timeStep *= 0.5;

        }
        else {
            /* Accept the step, increase the time step and update x*/

            currentTime += timeStep;
            timeStep *= 2;
            if (timeStep > 1) // Max time step
            {
                timeStep = 1;
            }

            // Update the old vectors depending on the order
            for (size_t i = 0; i < vectorOfNonLinearSystems.size(); i++)
            {
                for (int j = vectorOfNonLinearSystems[i]->order - 1; j > 0; j--)
                {
                    vectorOfNonLinearSystems[i]->xPrev[j] = vectorOfNonLinearSystems[i]->xPrev[j - 1];
                }
                vectorOfNonLinearSystems[i]->xPrev[0] = vectorOfNonLinearSystems[i]->x;
            }

            for (size_t i = 0; i < vectorOfNonLinearSystems.size(); i++)
            {
                VectorPrinter::printVector(vectorOfNonLinearSystems[i]->x);
            }

            std::cout << "error = " << error << "\tOutter Itterations = " << itterations << "\n";

        }

    }

};

#pragma once

#include "INonLinearSystem.h"

class ContinuityPressureSystem : public INonLinearSystem {
public:

    ContinuityPressureSystem(
        std::vector<double> _x,
        int _order) :
        INonLinearSystem(_x, _order) {}
    ~ContinuityPressureSystem() {}

    std::vector<double>* u = nullptr;

    void gaussSeidel(const double& dt) override {

    }

    double evaluateError(const double& dt) override {
        return 0;
    }

};

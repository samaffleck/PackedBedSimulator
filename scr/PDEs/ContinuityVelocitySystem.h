#pragma once

#include "INonLinearSystem.h"
#include "../SystemObjects/PackedBed.h"

class PackedBed;
class INonLinearSystem;

class ContinuityVelocitySystem : public INonLinearSystem {
public:

    ContinuityVelocitySystem(
        PackedBed* _bed,
        std::vector<double>& _x,
        int _order) :
        INonLinearSystem(_bed, _x, _order) {}
    ~ContinuityVelocitySystem() {}

    void updateRHS(const double& dt) override;
    void updateVelocity();
    void updateLinkCoefficients(const double& dt) override {};
    void correctVelocity();

private:

    const double alpha = 0.6;   // Damping factor

};

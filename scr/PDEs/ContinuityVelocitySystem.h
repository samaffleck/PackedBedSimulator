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

};

#pragma once

#include "INonLinearSystem.h"
#include "../SystemObjects/PackedBed.h"

class PackedBed;
class INonLinearSystem;

class AdvectionDiffusionSystem : public INonLinearSystem {
public:

    AdvectionDiffusionSystem(
        PackedBed* _bed,
        std::vector<double>& _x,
        int _order) :
        INonLinearSystem(_bed, _x, _order) {}
    ~AdvectionDiffusionSystem() {}

    void updateRHS(const double& dt) override;
    void updateLinkCoefficients(const double& dt) override {};

};

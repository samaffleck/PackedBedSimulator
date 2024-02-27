#pragma once

#include "INonLinearSystem.h"
#include "../SystemObjects/PackedBed.h"

class PackedBed;
class INonLinearSystem;

class ContinuityDensitySystem : public INonLinearSystem {
public:

    ContinuityDensitySystem(
        PackedBed* _bed,
        std::vector<double>& _x,
        int _order) :
        INonLinearSystem(_bed, _x, _order){}
    ~ContinuityDensitySystem() {}
    
    void updateRHS(const double& dt) override;

};

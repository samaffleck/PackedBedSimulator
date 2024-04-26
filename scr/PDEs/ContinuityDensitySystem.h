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
    void updateLinkCoefficients(const double& dt) override;
    void BDF_1(const double& dt);
    void BDF_1_Density(const double& dt);
    void BDF_2(const double& dt);
    void BDF_3(const double& dt);
    void BDF_4(const double& dt);
    void BDF_5(const double& dt);
    void BDF_6(const double& dt);
    void CN(const double& dt);
    void SIMPLE(const double& dt);
    void correctPressure();

};

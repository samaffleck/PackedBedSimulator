#pragma once

#include "IStep.h"
#include "../BoundaryConditions/BoundaryCondition_Constant.h"
#include "../BoundaryConditions/IBoundaryCondition.h"

class Pressurize : public IStep {
public:
	Pressurize(IBoundaryCondition* _inletPressure,
		IBoundaryCondition* _inletTemperature) : 
		inletPressure(_inletPressure),
		inletTemperature(_inletTemperature){}
	~Pressurize() {}

	IBoundaryCondition* inletPressure = nullptr;
	IBoundaryCondition* inletTemperature = nullptr;
	
	double inletDensityRHS(PackedBed* bed, const std::vector<double>& x, const std::vector<std::vector<double>>& xPrev, const double& dt) override;
	double outletDensityRHS(PackedBed* bed, const std::vector<double>& x, const std::vector<std::vector<double>>& xPrev, const double& dt) override;
	
	double inletVelocityRHS(PackedBed* bed) override;
	double outletVelocityRHS(PackedBed* bed) override;

	double inletTemperatureRHS(PackedBed* bed, const std::vector<double>& x, const std::vector<std::vector<double>>& xPrev, const double& dt) override;
	double outletTemperatureRHS(PackedBed* bed, const std::vector<double>& x, const std::vector<std::vector<double>>& xPrev, const double& dt) override;

	void updateBoundaryConditions(const double& currentTime) override;

};

#pragma once

#include "IStep.h"
#include "../BoundaryConditions/BoundaryCondition_Constant.h"
#include "../BoundaryConditions/IBoundaryCondition.h"

class FlowThrough : public IStep {
public:
	FlowThrough(IBoundaryCondition* _inletVelocity, 
		IBoundaryCondition* _outletPressure,
		IBoundaryCondition* _inletTemperature) :
		inletVelocity(_inletVelocity), 
		outletPressure(_outletPressure),
		inletTemperature(_inletTemperature){}
	~FlowThrough() {}

	IBoundaryCondition* outletPressure = nullptr;
	IBoundaryCondition* inletVelocity = nullptr;
	IBoundaryCondition* inletTemperature = nullptr;

	double inletDensityRHS(PackedBed* bed, const std::vector<double>& x, const std::vector<std::vector<double>>& xPrev, const double& dt) override;
	double outletDensityRHS(PackedBed* bed, const std::vector<double>& x, const std::vector<std::vector<double>>& xPrev, const double& dt) override;

	double inletVelocityRHS(PackedBed* bed) override;
	double outletVelocityRHS(PackedBed* bed) override;

	double inletTemperatureRHS(PackedBed* bed, const std::vector<double>& x, const std::vector<std::vector<double>>& xPrev, const double& dt) override;
	double outletTemperatureRHS(PackedBed* bed, const std::vector<double>& x, const std::vector<std::vector<double>>& xPrev, const double& dt) override;

	void updateBoundaryConditions(const double& currentTime) override;

};

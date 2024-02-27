#pragma once

#include "../SystemObjects/PackedBed.h"

#include <vector>

class PackedBed;

class IStep {
public:
	IStep() {}
	virtual ~IStep() {}

	const double R = 8.314;

	virtual double inletDensityRHS(PackedBed* bed, const std::vector<double>& x, const std::vector<std::vector<double>>& xPrev, const double& dt) = 0;
	virtual double outletDensityRHS() = 0;
	
	virtual double inletVelocityRHS(PackedBed* bed) = 0;
	virtual double outletVelocityRHS(PackedBed* bed) = 0;
	
	virtual double inletTemperatureRHS(PackedBed* bed, const std::vector<double>& x, const std::vector<std::vector<double>>& xPrev, const double& dt) = 0;
	virtual double outletTemperatureRHS(PackedBed* bed, const std::vector<double>& x, const std::vector<std::vector<double>>& xPrev, const double& dt) = 0;

	virtual void updateBoundaryConditions(const double& currentTime) = 0;

};

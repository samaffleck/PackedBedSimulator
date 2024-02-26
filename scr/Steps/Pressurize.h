#pragma once

#include "IStep.h"
#include "../BoundaryConditions/BoundaryCondition_Constant.h"
#include "../BoundaryConditions/IBoundaryCondition.h"

class Pressurize : public IStep {
public:
	Pressurize(IBoundaryCondition* _inletPressure)
	{
		inletPressure = _inletPressure;
	}
	~Pressurize() {}

	IBoundaryCondition* inletPressure = nullptr;
	
	double inletPressureRHS(const std::vector<std::vector<double>>& xPrev, const double& dt, const double& dx, const std::vector<double>& x, const std::vector<double>& u_x, const double& vis, const double& kappa) override;
	double outletPressureRHS() override;
	double inletVelocityRHS(const double& kappa, const double& dx, const double& vis, const std::vector<double>& p_x) override;
	double outletVelocityRHS(const double& kappa, const double& dx, const double& vis, const std::vector<double>& p_x, const int& lastIndex) override;

	void updateBoundaryConditions(const double& currentTime) override;

};

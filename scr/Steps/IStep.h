#pragma once

#include <vector>

class IStep {
public:
	IStep() {}
	virtual ~IStep() {}

	virtual double inletPressureRHS(const std::vector<std::vector<double>>& xPrev, const double& dt, const double& dx, const std::vector<double>& x, const std::vector<double>& u_x, const double& vis, const double& kappa) = 0;
	virtual double outletPressureRHS() = 0;
	virtual double inletVelocityRHS(const double& kappa, const double& dx, const double& vis, const std::vector<double>& p_x) = 0;
	virtual double outletVelocityRHS(const double& kappa, const double& dx, const double& vis, const std::vector<double>& p_x, const int& lastIndex) = 0;

	virtual void updateBoundaryConditions(const double& currentTime) = 0;

};

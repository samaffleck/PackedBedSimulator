#pragma once

#include "IBoundaryCondition.h"

class BoundaryCondition_Ramp :public IBoundaryCondition {
public:

	BoundaryCondition_Ramp(double _initialBoundaryValue, double _stepBoundaryValue, double _rampTime) :
		IBoundaryCondition(BoundaryConditionType::RAMP, _initialBoundaryValue),
		initialBoundaryValue(_initialBoundaryValue),
		stepBoundaryValue(_stepBoundaryValue),
		rampTime(_rampTime) {};
	~BoundaryCondition_Ramp() {};

	double initialBoundaryValue;
	double stepBoundaryValue;
	double rampTime;

	void update(const double& currentTime) override {
		if (currentTime < rampTime)
		{
			boundaryValue = initialBoundaryValue + ((stepBoundaryValue - initialBoundaryValue) / rampTime) * currentTime;
		}
		else {
			boundaryValue = stepBoundaryValue;
		}
	}

};
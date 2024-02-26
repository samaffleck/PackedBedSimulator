#pragma once

#include "IBoundaryCondition.h"

class BoundaryCondition_Step :public IBoundaryCondition {

public:

	BoundaryCondition_Step(double _initialBoundaryValue, double _stepBoundaryValue, double _startTime) :
		IBoundaryCondition(BoundaryConditionType::STEP, _initialBoundaryValue),
		initialBoundaryValue(_initialBoundaryValue),
		stepBoundaryValue(_stepBoundaryValue),
		startTime(_startTime){};
	~BoundaryCondition_Step() {};

	double initialBoundaryValue;
	double stepBoundaryValue;
	double startTime;

	void update(const double& currentTime) override {
		if (currentTime > startTime)
		{
			boundaryValue = stepBoundaryValue;
		}
	}

};
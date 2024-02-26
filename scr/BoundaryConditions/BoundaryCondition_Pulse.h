#pragma once

#include "IBoundaryCondition.h"

class BoundaryCondition_Pulse :public IBoundaryCondition {

public:

	BoundaryCondition_Pulse(double _initialBoundaryValue, double _stepBoundaryValue, double _startTime, double _endTime) :
		IBoundaryCondition(BoundaryConditionType::PULSE, _initialBoundaryValue),
		initialBoundaryValue(_initialBoundaryValue),
		stepBoundaryValue(_stepBoundaryValue),
		startTime(_startTime),
		endTime(_endTime){};
	~BoundaryCondition_Pulse() {};

	double initialBoundaryValue;
	double stepBoundaryValue;
	double startTime;
	double endTime;

	void update(const double& currentTime) override {
		if (currentTime > startTime and currentTime < endTime)
		{
			boundaryValue = stepBoundaryValue;
		}
		else {
			boundaryValue = initialBoundaryValue;
		}
	}

};


#pragma once

#include "IBoundaryCondition.h"

class BoundaryCondition_Constant :public IBoundaryCondition {

public:

	BoundaryCondition_Constant(double _boundaryValue) : 
		IBoundaryCondition(BoundaryConditionType::CONSTANT, _boundaryValue) {};
	~BoundaryCondition_Constant() {};

	void update(const double& currentTime) override {
		// No update;
	}

};

#pragma once

#include "IBoundaryCondition.h"

class BoundaryCondition_Flux :public IBoundaryCondition {

public:

	BoundaryCondition_Flux(double _boundaryValue) :
		IBoundaryCondition(BoundaryConditionType::FLUX, _boundaryValue) {};
	~BoundaryCondition_Flux() {};

	void update(const double& currentTime) override {
		// No update;
	}

};

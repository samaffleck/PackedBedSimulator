#pragma once

#include "IBoundaryCondition.h"

class BoundaryCondition_Flux :public IBoundaryCondition {

public:

	BoundaryCondition_Flux(double _boundaryValue) :
		IBoundaryCondition(BoundaryConditionType::FLUX, boundaryValue) {};
	~BoundaryCondition_Flux() {};

};

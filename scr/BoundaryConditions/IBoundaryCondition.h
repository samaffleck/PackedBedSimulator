#pragma once

enum class BoundaryConditionType {
	CONSTANT,
	FLUX
};

class IBoundaryCondition {
public:
	IBoundaryCondition(BoundaryConditionType _type, double _boundaryValue) : 
		type(_type),
		boundaryValue(_boundaryValue){}
	~IBoundaryCondition() {}

	BoundaryConditionType type;
	double boundaryValue;

};

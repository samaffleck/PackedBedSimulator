#pragma once

enum class BoundaryConditionType {
	CONSTANT,
	FLUX,
	STEP,
	PULSE,
	RAMP
};

class IBoundaryCondition {
public:
	IBoundaryCondition(BoundaryConditionType _type, double _boundaryValue) : 
		type(_type),
		boundaryValue(_boundaryValue){}
	virtual ~IBoundaryCondition() = default;

	BoundaryConditionType type;
	double boundaryValue;

	virtual void update(const double& currentTime) = 0;

};

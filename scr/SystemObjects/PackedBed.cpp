#include "PackedBed.h"


void PackedBed::initialise(const int& _numberOfCells, const double& initialTemperature, const double& initialVelocity, const double& initialPressure) {

	numberOfCells = _numberOfCells;
	dx = length / _numberOfCells;

	C.resize(_numberOfCells);
	U.resize(_numberOfCells);
	T.resize(_numberOfCells);

	for (int i = 0; i < _numberOfCells; i++)
	{
		C[i] = initialPressure / (R * initialTemperature);
		U[i] = initialVelocity;
		T[i] = initialTemperature;
	}

}


void PackedBed::selectStep(IStep* _step) {
	step = _step;
}

#include "PackedBed.h"


void PackedBed::initialise(const int& _numberOfCells, const double& initialTemperature, const double& initialVelocity, const double& initialPressure) {

	numberOfCells = _numberOfCells;
	dx = length / _numberOfCells;

	C.resize(_numberOfCells);
	U.resize(_numberOfCells + 1);
	T.resize(_numberOfCells);
	P.resize(_numberOfCells);
	P_old.resize(_numberOfCells);
	dP.resize(_numberOfCells);

	for (int i = 0; i < _numberOfCells; i++)
	{
		C[i] = (initialPressure + referencePressure) / (R * initialTemperature);
		U[i] = initialVelocity;
		T[i] = initialTemperature;
		P[i] = initialPressure;
		P_old[i] = initialPressure;
	}

	U[_numberOfCells] = initialVelocity;

}


void PackedBed::selectStep(IStep* _step) {
	step = _step;
}

void PackedBed::updateConstants()
{

	for (int n = 0; n < numberOfCells; n++)
	{
		P[n] = (C[n] * R * T[n]) - referencePressure;
		//C[n] = (P[n] + referencePressure) / (R * T[n]);
	}

}

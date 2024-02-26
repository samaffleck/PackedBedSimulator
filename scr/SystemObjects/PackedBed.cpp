#include "PackedBed.h"


void PackedBed::initialise(const int& _numberOfCells, const double& initialTemperature, const double& initialVelocity, const double& initialPressure) {

	numberOfCells = _numberOfCells;
	dx = length / _numberOfCells;

	P.resize(_numberOfCells);
	U.resize(_numberOfCells);
	T.resize(_numberOfCells);

	for (int i = 0; i < _numberOfCells; i++)
	{
		P[i] = initialPressure;
		U[i] = initialVelocity;
		T[i] = initialTemperature;
	}

	pressureSystem = ContinuityPressureSystem(kappa, viscosity, P, 1, _numberOfCells, length);
	velocitySystem = ContinuityVelocitySystem(kappa, viscosity, U, 1, _numberOfCells, length);

	pressureSystem.uSystem = &velocitySystem;
	velocitySystem.pSystem = &pressureSystem;

}

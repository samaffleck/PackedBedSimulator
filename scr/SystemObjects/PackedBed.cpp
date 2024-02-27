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

	densitySystem = ContinuityDensitySystem(kappa, viscosity, et, C, 1, _numberOfCells, length);
	velocitySystem = ContinuityVelocitySystem(kappa, viscosity, U, 1, _numberOfCells, length);
	temperatureSystem = AdvectionDiffusionSystem(T,1, _numberOfCells, length);

	densitySystem.uSystem = &velocitySystem;
	densitySystem.tSystem = &temperatureSystem;
	velocitySystem.cSystem = &densitySystem;
	velocitySystem.tSystem = &temperatureSystem;

}


void PackedBed::selectStep(IStep* step) {
	densitySystem.step = step;
	velocitySystem.step = step;
	temperatureSystem.step = step;
}

#pragma once

#include <vector>

#include "../PDEs/ContinuityPressureSystem.h"
#include "../PDEs/ContinuityVelocitySystem.h"
#include "../PDEs/DiffusionSystem.h"


class PackedBed {
public:

	PackedBed(
		double _length,
		double _bedVoidage,
		double _particleVoidage,
		double _particleDiameter,
		double _viscosity) :
		length(_length),
		eb(_bedVoidage),
		ep(_particleVoidage),
		dp(_particleDiameter),
		viscosity(_viscosity)
	{
		kappa = dp * dp * eb * eb * eb / (180 * (1 - eb) * (1 - eb));
		et = eb + (1 - eb) * ep;
	}
	~PackedBed() {}

	int numberOfCells{};
	double length;
	double eb;			// Bed void fraction
	double ep;			// Particle void fraction
	double et;			// Total voidage
	double dp;			// Particle diameter
	double kappa;		// Pereability
	double viscosity; 
	double order = 1;
	double dx{};

	// Initial conditions
	std::vector<double> P{};
	std::vector<double> U{};
	std::vector<double> T{};

	// Solver systems
	ContinuityPressureSystem pressureSystem;
	ContinuityVelocitySystem velocitySystem;

	void initialise(const int& systemSize, const double& initialTemperature, const double& initialVelocity, const double& initialPressure);

};

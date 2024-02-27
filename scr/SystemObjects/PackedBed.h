#pragma once

#include <vector>

#include "../Steps/IStep.h"

class IStep;

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
	int order = 1;
	double dx{};
	const double R = 8.314;
	double lambda = 1e-4;

	// Initial conditions
	std::vector<double> C{};	// Molar density [mol/m3]
	std::vector<double> U{};	// Velocity [m/s]
	std::vector<double> T{};	// Temperature [K]

	IStep* step = nullptr;

	void initialise(const int& systemSize, const double& initialTemperature, const double& initialVelocity, const double& initialPressure);
	void selectStep(IStep* _step);

};

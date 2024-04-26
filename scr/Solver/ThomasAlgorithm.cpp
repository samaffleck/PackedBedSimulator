#include <vector>
#include <iostream>

#include "ThomasAlgoritms.h"

void ThomasAlgorithm::initialise(int size) {
	n = size;
	cPrime.resize(n);
	dPrime.resize(n);
}


void ThomasAlgorithm::solve(const std::vector<double>& Ae, const std::vector<double>& Ao, const std::vector<double>& Aw, const std::vector<double>& X, std::vector<double>& V) {
	// AV = X

	// Forward pass
	cPrime[0] = Ae[0] / Ao[0];
	dPrime[0] = X[0] / Ao[0];
	for (int i = 1; i < n - 1; i++)
	{
		cPrime[i] = Ae[i] / (Ao[i] - Aw[i] * cPrime[i - 1]);
		dPrime[i] = (X[i] - Aw[i] * dPrime[i - 1]) / (Ao[i] - Aw[i] * cPrime[i - 1]);
	}
	dPrime[n - 1] = (X[n - 1] - Aw[n - 1] * dPrime[n - 2]) / (Ao[n - 1] - Aw[n - 1] * cPrime[n - 2]);
	
	// Back substitution
	V[n - 1] = dPrime[n - 1];
	for (int i = n - 2; i >= 0; i--)
	{
		V[i] = dPrime[i] - cPrime[i] * V[i + 1];
	}

}


void ThomasAlgorithm::printVector(const std::vector<double>& v) {
	std::cout << "[" << v[0];
	for (size_t i = 1; i < v.size(); i++)
	{
		std::cout << ", " << v[i];
	}
	std::cout << "]";
}

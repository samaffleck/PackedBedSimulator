#pragma once

#include <vector>

class ThomasAlgorithm {

public:

	ThomasAlgorithm() {};
	~ThomasAlgorithm() {};

	int n{};
	std::vector<double> cPrime{};
	std::vector<double> dPrime{};

	void initialise(int size);
	void solve(const std::vector<double>& Ae, const std::vector<double>& Ao, const std::vector<double>& Aw, const std::vector<double>& X, std::vector<double>& V);
	void printVector(const std::vector<double>& v);

};

#pragma once

#include <fstream>
#include <vector>

#include "../PDEs/INonLinearSystem.h"


class DataLogger {
public:
	DataLogger() {}
	~DataLogger() {}

	std::vector<std::ofstream> vectorOfFiles{};

	void open(const int& numberOfFiles, const int& numberOfNodes, const double& dx);
	void close();
	void log(const int& systemNumber, const double& time, const std::vector<double>& data);

};

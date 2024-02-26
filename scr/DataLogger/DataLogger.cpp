#include "DataLogger.h"
#include <filesystem>
#include <string>

void DataLogger::open(const int& numberOfFiles, const int& numberOfNodes, const double& dx) {

    // Check if results folder exists and create if not
    std::filesystem::create_directory("results");

    for (int i = 0; i < numberOfFiles; i++)
    {
        vectorOfFiles.push_back(std::ofstream());
        vectorOfFiles[i].open("results/system_" + std::to_string(i) + ".csv");
        vectorOfFiles[i] << "t [s] | x [m]]";

        for (int n = 0; n < numberOfNodes; n++)
        {
            vectorOfFiles[i] << "," << n * (dx);
        }

        vectorOfFiles[i] << std::endl;

    }

}

void DataLogger::close() {

    for (size_t i = 0; i < vectorOfFiles.size(); i++)
    {
        vectorOfFiles[i].close();
    }

}

void DataLogger::log(const int& systemNumber, const double& time, const std::vector<double>& data) {

    vectorOfFiles[systemNumber] << time;

    for (size_t n = 0; n < data.size(); n++)
    {
        vectorOfFiles[systemNumber] << "," << data[n];
    }

    vectorOfFiles[systemNumber] << std::endl;

}

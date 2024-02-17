#pragma once

#include <vector>


class VectorNorm {
public:

    static double L2Norm(const std::vector<double> e) {
        double error = 0;
        for (auto itt = e.begin(); itt != e.end(); itt++)
        {
            error += (*itt) * (*itt);
        }

        return sqrt(error);
    }

private:
    VectorNorm() {}
    ~VectorNorm() {}

};

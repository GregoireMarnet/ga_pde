#ifndef VOLATILITY_HPP
#define VOLATILITY_HPP

#include <iostream>
#include <vector>
#include "eigen-3.4.0/Eigen/Core"

namespace dauphine {

class volatility{

    public : 
    //Eigen::MatrixXd get_vol();
    //double operator() (const int&t, const int&x);
    //std::vector<double> operator() (const int&t);

    protected:
    //Eigen::MatrixXd surface;
};

class constant_vol : public volatility
{
public:
    constant_vol(double sigma);
};

};

#endif
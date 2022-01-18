#ifndef VOLATILITY_HPP
#define VOLATILITY_HPP

#include <iostream>
#include <vector>
#include "mesh.hpp"


namespace dauphine {

class volatility{

    public : 
        volatility();
        std::vector<double> get_vol() const;

    protected:
        std::vector<double> p_vol;
};

class vol_BS : public volatility
{
    public:
        vol_BS(double sigma);
};

class vol_gen : public volatility
{
    public:
        vol_gen(double sigma, dauphine::mesh msh,  double alpha);
    
};

};

#endif
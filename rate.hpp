#ifndef rate_hpp
#define rate_hpp

#include <vector>
#include "mesh.hpp"

namespace dauphine

{

class rate 
{
    public:
        rate();
        std::vector<double> get_rates() const;

    protected:
        std::vector<double> p_rate;
};

class rate_BS : public rate
{
    public:
        rate_BS(double r);
};

class rate_gen : public rate
{
    public:
        rate_gen(double r, dauphine::mesh msh,  double alpha);
    
};

}

#endif
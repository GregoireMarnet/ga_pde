#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP

#include "payoff.hpp"

namespace dauphine

{
class boundary
{
    public:
        boundary();
        virtual double get_upper_b() const =0;
        virtual double get_lower_b() const =0;

    protected:
        double upper_b;
        double lower_b;        
};

class dirichlet : public boundary
{
    
    public:
        dirichlet(dauphine::payoff payoff, double s_min, double s_max);
        virtual double get_upper_b() const;
        virtual double get_lower_b() const;
    
};

}


#endif
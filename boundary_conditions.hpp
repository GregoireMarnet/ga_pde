#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP

namespace dauphine

{
class boundary
{
    public:
        boundary();
        double get_upper_b() const =0;
        double get_lower_b() const =0;

    protected:
        double upper_b;
        double lower_b;        
};

class dirichlet : public boundary
{
    
    public:
        dirichlet(payoff::payoff payoff, double s_min, double s_max);
        double get_upper_b() const;
        double get_lower_b() const;
    
};

}


#endif
#include "boundary_conditions.hpp"


namespace dauphine

{

boundary::boundary(){}

dirichlet::dirichlet(dauphine::payoff payoff, double s_min, double s_max)
{
    upper_b = payoff(s_max);
    lower_b = payoff(s_min);
}

double get_upper_b() const
{
    return upper_b;
}

double get_lower_b() const
{
    return lower_b;
}

}
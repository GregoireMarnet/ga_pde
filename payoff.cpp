#include "payoff.hpp"
#include <algorithm> //for std::max function


namespace payoff
{
payoff::payoff(){} 

payoff_call::payoff_call(const double& strike)
{
    p_strike = strike;
}

double payoff_call::operator() (const double& spot) const
{
    return std::max(spot - p_strike, 0.0);
}


}
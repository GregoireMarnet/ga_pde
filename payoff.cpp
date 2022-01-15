#include "payoff.hpp"
#include <algorithm> //for std::max function
#include <math.h>

namespace dauphine

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

std::vector<double> payoff_call::operator() ( const std::vector<double>& spot) const
{
    std::vector<double> payout(spot.size());
    for( int i=0; i<spot.size();i++){
        payout[i]=std::max(exp(spot[i]) - p_strike, 0.0);
    }
    return payout;
}

}
#include "payoff.hpp"


payoff::payoff(const std::function<double(double)>& function_spot)
    : p_payoff(function_spot)
{
}

std::function<double(double)> function()
{
    return p_payoff;
}
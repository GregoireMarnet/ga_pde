#include "payoff.hpp"


payoff::payoff(const std::function<double(double)>& function_spot)
    : p_payoff(function_spot)
{
}
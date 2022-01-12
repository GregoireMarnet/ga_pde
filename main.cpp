#include "closed_form.hpp"
#include "payoff.hpp"
#include "boundary_conditions.hpp"
#include <iostream>


// Guidelines:
//
// 1] Equation
// Even if the volatility and the rate are constant in the BS model,
// we can prove that, under certain assumptions, we can use the same
// PDE equation with volatility surface and rate curves. We could also
// want to take into account the repo (even if it could theoretically
// be part of the r factor). Therefore, it is more generic to solve
// the following generic equation;
// df/dt = a(x,t)d2f/dx2 + b(x,t)df/dx + c(x, t)f + d(x, t).
// The pricer should ask the coefficients a, b, c and d to an
// abstract class that will be inherited by classes implementing
// different models.
// 
// 2] Payoff
// The pricer should be able to price exotic options, such as
// barriers or asian options. To do so, the pricer must be able
// to call an external function between each step. Define an API
// that allows to register a function or an abstract class modeling
// a payoff.

//# include <eigen3/Eigen>
#include "eigen-3.4.0/Eigen/Core"

int main(int argc, const char * argv[])
{
    const double strike = 20;
    const double spot = 23;

    dauphine::payoff_call payoff_call(strike);    
    
    std::cout << "Strike is : " << strike << std::endl;
    std::cout << "Spot is : " << spot << std::endl;
    std::cout << "Payoff is :" << payoff_call(spot) << std::endl;  

    dauphine::dirichlet bound(payoff_call,0,30);

    std::cout << "Upper Bound : " << bound.get_upper_b() << std::endl;
    std::cout << "Lower Bound : " << bound.get_lower_b() << std::endl;


    
    return 0;
}

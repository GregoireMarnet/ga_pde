#include "closed_form.hpp"
#include "payoff.hpp"
#include "boundary_conditions.hpp"
#include "volatility.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "rate.hpp"


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

std::ostream &operator<<(std::ostream &os, const std::vector<double> &input)
{
    for (auto const &i: input) {
        os << i << " ";
    }
    return os;
}

int main(int argc, const char * argv[])
{
    const double strike = 20;
    const double spot = 23;

    dauphine::payoff_call poff_call(strike);    
    
    std::cout << "Strike is : " << strike << std::endl;
    std::cout << "Spot is : " << spot << std::endl;
    std::cout << "Payoff is :" << poff_call(spot) << std::endl;  

    dauphine::dirichlet bound(poff_call,0,30);

    std::cout << "Upper Bound : " << bound.get_upper_b() << std::endl;
    std::cout << "Lower Bound : " << bound.get_lower_b() << std::endl;

    const double maturity = 10;

    dauphine::mesh msh(spot, maturity,3,10,0.05);
    /*dauphine::solver solv;

    std::cout << "xaxis : " << msh.get_xaxis() << std::endl;
    std::cout << "final_cond: " << solv.compute_price(poff_call,msh, bound) << std::endl;*/



    const double r0 = 3.;
    const double alpha = 0.1;

    dauphine::rate_BS test_rate1(r0);
    dauphine::rate_gen test_rate2(r0,msh, alpha);

    std::cout << "Rate : " << test_rate1.get_rates() << "%" << std::endl;
    std::cout << "Rate : " << test_rate2.get_rates() << std::endl;





    return 0;
}

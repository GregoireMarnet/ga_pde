#include "closed_form.hpp"
#include "payoff.hpp"
#include "boundary_conditions.hpp"
#include "volatility.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "rate.hpp"
#include "matrix.hpp"


#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h> 
//#include "eigen-3.4.0/Eigen/Core"


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



std::ostream &operator<<(std::ostream &out, const std::vector<double> &input)
{
    for(std::size_t i = 0; i < input.size(); ++i)
        {
                out << input[i] << " ";
        }
        out << std::endl;
        return out;
    return out;
}

int main(int argc, const char * argv[])
{
    /*
    const double strike = 20;
    const double spot = 23;

    dauphine::payoff_call poff_call(strike);    
    
    std::cout << "Strike is : " << strike << std::endl;
    std::cout << "Spot is : " << spot << std::endl;
    std::cout << "Payoff is :" << poff_call(spot) << std::endl;  

    dauphine::dirichlet bound(poff_call,0,30);

    std::cout << "Upper Bound : " << bound.get_upper_b() << std::endl;
    std::cout << "Lower Bound : " << bound.get_lower_b() << std::endl;

    const double maturity = 1;
    const double r0 = 0.03;
    const double v0 = 0.2;


    dauphine::vol_BS vol(v0);
    dauphine::rate_BS rate(r0);
    dauphine::mesh msh(spot, maturity,11,3,v0);

    dauphine::solver solv(poff_call,msh, bound,vol, rate, 0.5);
    solv.compute_price();

    */


    std::vector<double> vect = {1,2,3,4,5};
    std::vector<double> vect2 = {2,2,3,4,5};
    //std::copy(vect2.begin(),vect2.end(),vect.begin()+1);

    //std::transform(vect.begin(),vect.end(),vect.begin(),[](double f)-> double {return f*f;});

    std::cout << vect - vect2;
    return 0;
}

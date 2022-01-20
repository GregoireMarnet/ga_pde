#include "closed_form.hpp"
#include "payoff.hpp"
#include "boundary_conditions.hpp"
#include "volatility.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "rate.hpp"
#include "matrix.hpp"
#include "pde.hpp"


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
    
    const double strike = 100;
    const double spot = 100;
    const double maturity = 1;
    const double r0 = 0.05;
    const double v0 = 0.1;
    const int ndx = 1001;
    const int ndt = 1500;
    const double theta = 0.5;

    const double kappa = 2;
    const double heston_theta = 0.12;

    
    dauphine::vol_BS vol(v0);
    dauphine::rate_BS rate(r0);
    dauphine::mesh msh(spot, maturity,ndx,ndt,v0);
    dauphine::payoff_call poff_call(strike); 
    dauphine::dirichlet bound(poff_call,msh.get_xmin(),msh.get_xmax());
    dauphine::pde_european_BS pde_BS(vol,rate);
    

    dauphine::vol_heston vol_h(v0,kappa,heston_theta,msh);
    dauphine::pde_european_heston pde_h(vol_h,rate);

    dauphine::solver solv(poff_call,msh, bound, theta);
    
    std::cout << solv << std::endl;
    solv.call_compute_price(pde_BS);
    solv.call_compute_price(pde_h);



    





    
    
    return 0;
}

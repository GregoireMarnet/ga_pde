#include "closed_form.hpp"
#include "payoff.hpp"
#include "boundary_conditions.hpp"
#include "volatility.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "rate.hpp"
#include "matrix.hpp"
#include "pde.hpp"
#include "greeks.hpp"
#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h> 


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

 /* print vector */
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
    // Parameters of the option
    const bool is_call = true;
    const double spot = 100.;
    const double strike = 100.;
    const double maturity = 1.; // in years
    const double v0 = 0.10;
    const double r0 = 0.05;
    
    // Parameters of the mesh space 
    const int ndx = 501; // enter odd number 
    const int ndt = 100;

    // Theta scheme parameter
    const double theta = 0.5;

    // Parameters to implement an heston model 
    const double kappa = 2;
    const double heston_theta = 0.12;

    // Mesh space definition
    dauphine::mesh msh(spot, maturity,ndx,ndt,v0);
    dauphine::vanilla_poff poff_call(strike, is_call); 
    dauphine::dirichlet bound(poff_call,msh.get_xmin(),msh.get_xmax());

    // Constant parameters definition
    dauphine::vol_BS vol(v0);
    dauphine::rate rate(r0,msh,0.);

    // Non-constant parameters definition
    dauphine::vol_heston vol_h(v0,kappa,heston_theta,msh);
    dauphine::rate rg(r0, msh, 0.);
    
    // Definition of the pde to solve
    dauphine::pde_european pde(vol,rate);
    
    // Building finite difference method parameters
    dauphine::solver solv(poff_call, msh, bound, theta);
    std::cout << solv << std::endl;

    // Results
    double cf  =  dauphine::bs_price(spot,strike,v0,maturity,is_call) ;
    std::cout << "Closed form price : " << cf << " (with assumption r=0)" << std::endl;
    dauphine::matrix result = solv.price(pde);
    std::cout << "Finite difference method price : " << result((ndx-1)/2,0) << std::endl << std::endl;
    
    //Greeks
    double delta = dauphine::delta(result,msh.get_xaxis(),ndx/2);
    std::cout << "Delta : " << delta << std::endl;

    double gamma = dauphine::gamma(result,msh.get_xaxis());
    std::cout << "Gamma : " << gamma << std::endl;

    double rho = dauphine::rho(result,solv,vol,rate, ndx);
    std::cout << "Rho : " <<  rho << std::endl;

    double vega = dauphine::vega(result,solv,vol,rate, ndx);
    std::cout << "Vega : " <<  vega << std::endl;

    double t = dauphine::theta(result,msh.get_dt(),ndx/2);
    std::cout << "Theta : " << t << std::endl;

    return 0;
}

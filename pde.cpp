#include "pde.hpp"

namespace dauphine 
{
    pde::pde(){};

    pde_european::pde_european(){};

    double pde_european::get_coeff_a(){
        return p_a;
    }

    double pde_european::get_coeff_b(){
        return p_b;
    }

    double pde_european::get_coeff_c(){
        return p_c;
    }

    double pde_european::get_coeff_d(){
        return p_d;
    }
    
    pde_european_BS::pde_european_BS(const vol_BS& vol,const rate_BS& rate)
    {
        p_a = - 0.5 * pow(vol.get_vol()[0],2);
        p_b = 0.5 * pow(vol.get_vol()[0],2) - rate.get_rates()[0];
        p_c = rate.get_rates()[0];
        p_d = 0;
    }
}
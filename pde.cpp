#include "pde.hpp"

namespace dauphine 
{
    pde::pde(){};

    pde_european::pde_european(){};

    std::vector<double> pde_european::get_coeff_a() const{
        return p_a;
    }

    std::vector<double> pde_european::get_coeff_b() const{
        return p_b;
    }

    std::vector<double> pde_european::get_coeff_c() const{
        return p_c;
    }

    std::vector<double> pde_european::get_coeff_d() const{
        return p_d;
    }
    
    pde_european_BS::pde_european_BS(const vol_BS& vol,const rate_BS& rate)
    {
        p_a.resize(1);
        p_b.resize(1);
        p_c.resize(1);
        p_d.resize(1);
        p_a[0] = - 0.5 * pow(vol.get_vol()[0],2);
        p_b[0] = 0.5 * pow(vol.get_vol()[0],2) - rate.get_rates()[0];
        p_c[0] = rate.get_rates()[0];
        p_d[0]=0;
    }
    
    pde_european_gen::pde_european_gen(const vol_gen& vol, const rate_gen& rate)
    {
        const int size = vol.get_vol().size();
        p_a.resize(size);
        p_b.resize(size);
        p_c.resize(size);
        p_d.resize(size);
        for (int i=0; i<size; i++)
        {
            p_a[i] = - 0.5 * pow(vol.get_vol()[i],2);
            p_b[i] = 0.5 * pow(vol.get_vol()[i],2) - rate.get_rates()[i];
            p_c[i] = rate.get_rates()[i];
            p_d[i]=0; 
        }
    }

}
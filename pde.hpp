#ifndef PDE_HPP
#define PDE_HPP

#include "volatility.hpp"
#include "rate.hpp"
#include <math.h>

namespace dauphine 
{
    class pde
    {
        public:
            pde();
    };

    class pde_european : public pde // virtuelle pure
    {
        public :
            pde_european();
            double get_coeff_a();
            double get_coeff_b();
            double get_coeff_c();
            double get_coeff_d();

        protected : 
            double p_a;
            double p_b;
            double p_c;
            double p_d;

    };

        class pde_european_BS : public pde_european
        {
            public:
                pde_european_BS(const vol_BS& vol, const rate_BS& rate);
        };





}


#endif
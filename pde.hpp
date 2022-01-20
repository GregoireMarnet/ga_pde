#ifndef PDE_HPP
#define PDE_HPP

#include "volatility.hpp"
#include "rate.hpp"
#include <math.h>
#include <algorithm>

namespace dauphine 
{
    class pde
    {
        public:
            pde();
    };

    class pde_european : public pde 
    {
        public :
            pde_european(const volatility& vol, const rate& rate);
            std::vector<double> get_coeff_a() const;
            std::vector<double> get_coeff_b() const;
            std::vector<double> get_coeff_c() const;
            std::vector<double> get_coeff_d() const;
        protected : 
            std::vector<double> p_a;
            std::vector<double> p_b;
            std::vector<double> p_c;
            std::vector<double> p_d;

    };

        class pde_european_BS : public pde_european
        {
            public:
                pde_european_BS(const vol_BS& vol, const rate_BS& rate);
                
            
        };

        
        class pde_european_gen : public pde_european{
            public:
                pde_european_gen(const volatility& vol, const rate& rate);
        };

        class pde_european_heston : public pde_european{
            public:
                pde_european_heston(const vol_heston& vol, const rate_BS& rate);
        };
        






}


#endif
#include "volatility.hpp"

namespace dauphine
{

volatility::volatility(){};

vol_BS::vol_BS(double sigma)
{
    p_vol.resize(1);
    p_vol[0] = sigma;
}

vol_gen::vol_gen(double sigma, dauphine::mesh msh,  double alpha)
{
    const int n_dt = msh.get_ndt();
    p_vol.resize(n_dt);
    p_vol[0] = sigma;
    for(int i=1; i<n_dt ; i++){
            p_vol[i] = p_vol[i-1] + alpha;
        };
}

}
#include "rate.hpp"


namespace dauphine
{

rate::rate(){};

rate_BS::rate_BS (double r)
{   
    p_rate.resize(1);
    p_rate[0] = r;
}


rate_gen::rate_gen (double r, dauphine::mesh msh, double alpha)
{
    const int n_dt = msh.get_ndt();
    p_rate.resize(n_dt);
    p_rate[0] = r;
    for(int i=1; i<n_dt ; i++){
            p_rate[i] = p_rate[i-1] + alpha;
        };
}


std::vector<double> rate::get_rates() const
{
    return p_rate;
}


}
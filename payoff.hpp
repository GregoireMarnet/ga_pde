#ifndef PAYOFF_HPP
#define PAYOFF_HPP

namespace payoff
{

class payoff
{
    public :

        payoff(const std::function<double(double)>& function_spot);

    protected :

        std::function<double(double)> p_payoff;
};


};



#endif
#ifndef PAYOFF_HPP
#define PAYOFF_HPP

namespace payoff
{

double function_test(double spot)
    {
        return spot*spot;
    };

class payoff
{
    public :

        payoff(const std::function<double(double)>& function_spot);
        std::function<double(double)> function();

    protected :

        std::function<double(double)> p_payoff;
};


};



#endif
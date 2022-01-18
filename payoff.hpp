#ifndef PAYOFF_HPP
#define PAYOFF_HPP

#include <vector>

namespace dauphine

{

class payoff
{
    public :
        payoff();
        virtual double operator() (const double& spot) const = 0;
        virtual std::vector<double> operator() (const std::vector<double>& spot) const = 0;
        virtual double get_strike() const = 0;
};

class payoff_call : public payoff
{
    private : 
        double p_strike;

    public:
        payoff_call(const double& strike);
        virtual double get_strike() const;
        virtual double operator() (const double& spot) const;
        virtual std::vector<double> operator() (const std::vector<double>& spot) const;

};

};



#endif
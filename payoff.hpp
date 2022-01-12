#ifndef PAYOFF_HPP
#define PAYOFF_HPP


namespace dauphine

{

class payoff
{
    public :
        payoff();
        virtual double operator() (const double& spot) const = 0;
};

class payoff_call : public payoff
{
    private : 
        double p_strike;

    public:
        payoff_call(const double& strike);
        virtual double operator() (const double& spot) const;

};

};



#endif
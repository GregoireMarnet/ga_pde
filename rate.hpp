#ifndef rate_hpp
#define rate_hpp

#include <vector>


namespace dauphine

{

class rate 
{
    public:
        rate();
        double operator() (const int& i);

    protected:
        std::vector<double> rate_vector;
}

class rate_BS : public rate
{
    public:
        rate_BS(double r);

}

}

#endif
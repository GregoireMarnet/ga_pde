#include "rate.hpp"


namespace dauphine
{

rate::rate(){};

double rate::operator() (int i)
{
    return rate_vector[i]
}


}
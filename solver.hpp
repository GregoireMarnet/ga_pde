#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <vector>
#include "eigen-3.4.0/Eigen/Core"


namespace solver
{
std::vector<double> solver_mesh(payoff& payoff, 
                                const Eigen::Matrix<double,2,2>& a,
                                const std::vector<double>& b,
                                const std::vector<double>& c,
                                const std::vector<double>& d,
                                const double& theta,
                                const double& T,
                                const double& dt,
                                const double& s_min,
                                const double& s_max,
                                const double& dx,
                                dirichlet& bd);



void print_solver_mesh(double rate, double maturity);

}

#endif
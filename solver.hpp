#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <vector>
#include "payoff.hpp"
#include "eigen-3.4.0/Eigen/Core"


namespace solver
{

std::vector<double> solver_mesh(payoff::payoff& payoff, 
                                const Eigen::MatrixXd& a,
                                const Eigen::MatrixXd& b,
                                const Eigen::MatrixXd& c,
                                const Eigen::MatrixXd& d,
                                const double& theta,
                                const double& T,
                                const int& n_dt,
                                const double& s_min,
                                const double& s_max,
                                const int& n_dx,
                                boundary_condition& bd);



void print_solver_mesh(double rate, double maturity);

}

#endif
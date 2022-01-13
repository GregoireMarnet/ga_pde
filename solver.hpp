#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <vector>
#include "payoff.hpp"
#include "mesh.hpp"
#include "volatility.hpp"
#include "boundary_conditions.hpp"
#include "eigen-3.4.0/Eigen/Core"


namespace dauphine
{

    class solver {
        public : 
        std::vector<double> compute_price(payoff& poff, mesh &msh, boundary& bd);

        private : 
        /*mesh& m_mesh;
        volatility& m_v;
        payoff& m_p;
        boundary& m_bd;*/
    };

    std::vector<double> compute_price(payoff& poff, mesh &msh, boundary& bd);

/*std::vector<double> solver_mesh(dauphine::payoff& payoff, 
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
                                dauphine::boundary& bd);




void print_solver_mesh(double rate, double maturity);
*/
}

#endif
#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <vector>
#include "payoff.hpp"
#include "mesh.hpp"
#include "volatility.hpp"
#include "boundary_conditions.hpp"
#include "rate.hpp"
#include "eigen-3.4.0/Eigen/Core"


namespace dauphine
{

    class solver {
        public : 

        solver(payoff& poff,
                mesh& msh, 
                boundary& bd,
                volatility& vol,
                rate& rate,
                double theta);


        std::vector<double> compute_price();

        
        void init_coeff(double& a, double& b, double& c, double& d);
        void init_coeff(std::vector<double>& a,
                        std::vector<double>& b,
                        std::vector<double>& c,
                        std::vector<double>& d);

        void transform_coeff(double& a, double& b, double& c, double& d);
        void transform_coeff(std::vector<double>& a,
                        std::vector<double>& b,
                        std::vector<double>& c,
                        std::vector<double>& d);

        void init_matrice_1(Eigen::MatrixXd& m_trans, const int& dim, double a, double b, double c, double d);
        void init_matrice_2(Eigen::MatrixXd& m_trans, const int& dim, double a, double b, double c, double d);

        private :
            mesh& m_msh;
            volatility& m_vol;
            payoff& m_poff;
            boundary& m_bd;
            rate& m_rate;
            double m_theta;

    };

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
#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "payoff.hpp"
#include "mesh.hpp"
#include "volatility.hpp"
#include "boundary_conditions.hpp"
#include "rate.hpp"
#include "matrix.hpp"
#include <vector>
#include <ostream>
#include <iostream>
#include <math.h>


namespace dauphine
{

    class solver {
        public : 

        friend std::ostream& operator<<(std::ostream& os, solver const & s) {
            

            //os << "Strike is : " << s.m_poff.get_strike() << std::endl;
            os << "Spot is : " << s.m_msh.get_spot() << std::endl;
            //os << "Payoff is :" << s.m_poff(s.m_msh.get_spot()) << std::endl ;  
            os << std::endl;

            os << "Log(Smax)  : " << s.m_msh.get_xmax() << std::endl;
            os << "Log(Smin)  : " << s.m_msh.get_xmin() << std::endl;
            os << "Smax  : " << exp(s.m_msh.get_xmax()) << std::endl;
            os << "Smin  : " << exp(s.m_msh.get_xmin()) << std::endl;
            os << std::endl;

            os << "Upper Bound : " << s.m_bd.get_upper_b() << std::endl;
            os << "Lower Bound : " << s.m_bd.get_lower_b() << std::endl << std::endl;
            os << std::endl;

            os << "dt : " << s.m_msh.get_dt() << std::endl;
            os << "dx : " << s.m_msh.get_dx() << std::endl;
            os << std::endl;
    
        return os ;
         
        }

        solver(payoff& poff, mesh& msh, boundary& bd,volatility& vol,rate& rate,double theta);

        void compute_price();

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

        void init_matrice_1(matrix& m_trans, const int& dim, double a, double b, double c, double d);
        void init_matrice_2(matrix& m_trans, const int& dim, double a, double b, double c, double d);
        void fill_matrix(matrix& mesh_matrix, int t, std::vector<double> vect);

        std::vector<double> solve_tridiag(const matrix& mat, std::vector<double>& d);

        
        /*std::vector<double> solve_system(dauphine::matrix m1_inv, dauphine::matrix m_trans_2, 
                                        std::vector<double> final_vect,
                                        std::vector<double> vect, int ndx);

        */

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
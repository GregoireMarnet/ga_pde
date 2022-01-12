#include "solver.hpp"


namespace solver
{
    
std::vector<double> solver_mesh(payoff& payoff, 
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
                                boundary_condition& bd) // need to define it

{
    Eigen::Matrix<double, n_dx + 1, n_dt + 1> A;
    Eigen::Matrix<double, n_dx + 1, n_dt + 1> B;
    Eigen::Matrix<double, n_dx + 1, n_dt + 1> C;
    Eigen::Matrix<double, n_dx + 1, n_dt + 1> D;

    const double dt = T / n_dt;
    const double dx = (s_max - s_min) / n_dx;

    const double nu1 = dt / dS^2;
    const double nu2 = dt / dS;

    for ( int i=0 ; i < n_dx + 1 ; i++) {
	    for ( int j=0 ; j < n_dt + 1 ; j++ {
		    A(i,j) = theta * (a(i,j)*nu1 - (1/2) * b(i,j) * nu2) ;
		    B(i,j) = theta * (c(i,j) * dt - 2 * a(i,j) * nu1) ;
		    C(i,j) = theta * (a(i,j) * nu1  + (1/2) * b(i,j) * nu2) ;
		    D(i,j) = theta * d(i,j) * dt ;
	    }
    }


};

}

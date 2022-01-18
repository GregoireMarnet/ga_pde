#include "solver.hpp"
#include <math.h>
#include <algorithm>



namespace dauphine
{

    std::ostream &operator<<(std::ostream &out, const std::vector<double> &input)
    {
        for(std::size_t i = 0; i < input.size(); ++i)
            {
                    out << input[i] << " ";
            }
            out << std::endl;
            return out;
        return out;
    }

    solver::solver(payoff& poff,
                    mesh& msh, 
                    boundary& bd,
                    volatility& vol,
                    rate& rate,
                    double theta)
    : m_poff(poff), m_msh(msh), m_bd(bd), m_vol(vol), m_rate(rate), m_theta(theta)
    {};


    void solver::init_coeff(std::vector<double>& a,
                        std::vector<double>& b,
                        std::vector<double>& c,
                        std::vector<double>& d)
    {
        /*std::vector<double> var = m_vol.get_vol();
        std::transform(var.begin(),var.end(),var.begin(),[](double f)-> double {return f*f;});
        a = - 0.5 * var;
        b = 0.5 * var - m_rate.get_rates();
        c = m_rate.get_rates();
        d = 0;*/
    };

    void solver::init_coeff(double& a, double& b, double& c, double& d)
    {
        a = - 0.5 * pow(m_vol.get_vol()[0],2);
        b = 0.5 * pow(m_vol.get_vol()[0],2) - m_rate.get_rates()[0];
        c = m_rate.get_rates()[0];
        d = 0;
    };

    void solver::transform_coeff(std::vector<double>& a,
                        std::vector<double>& b,
                        std::vector<double>& c,
                        std::vector<double>& d)
    {
    };

    void solver::transform_coeff(double& a, double& b, double& c, double& d)
    {
        const double nu1 = m_msh.get_dt() / pow(m_msh.get_dx(),2);
        const double nu2 = m_msh.get_dt() / m_msh.get_dx();
        //std::cout << "nu1 : " << nu1 << " and nu2 : " << nu2 << std::endl;

        double A = a;
        double B = b;
        double C = c;
        double D = d;

        a = A*nu1 - 0.5 * B * nu2;
        b = C * m_msh.get_dt() - 2 * A * nu1;
        c = A * nu1  + 0.5 * B * nu2;
        d = D * m_msh.get_dt();
    };

    void solver::init_matrice_1(matrix& m_trans,const int& dim, double a, double b, double c, double d)
    {
        m_trans(0,0) = m_theta*b+1;
        m_trans(0,1) = m_theta*c;

        for(int i=1; i<dim-1; i++){
            m_trans(i,i-1) = m_theta*a;
            m_trans(i,i) = m_theta*b+1;
            m_trans(i,i+1) = m_theta*c;
        };

        m_trans(dim-1,dim-2) = m_theta*a;
        m_trans(dim-1,dim-1) = m_theta*b+1;

        
    }

    void solver::init_matrice_2(matrix& m_trans,const int& dim, double a, double b, double c, double d)
    {
        m_trans(0,0) = 1-(1-m_theta)*b;
        m_trans(0,1) = -(1-m_theta)*c;

        for(int i=1; i<dim-1; i++){
            m_trans(i,i-1) = -(1-m_theta)*a;
            m_trans(i,i) = 1-(1-m_theta)*b;
            m_trans(i,i+1) = -(1-m_theta)*c;
        };

        m_trans(dim-1,dim-2) = -(1-m_theta)*a;
        m_trans(dim-1,dim-1) = 1-(1-m_theta)*b;        
    }

    void solver::fill_matrix(matrix& mesh_matrix, int t, std::vector<double> vect)
    {
        for (int x=0; x<mesh_matrix.nb_rows();x++){
            mesh_matrix(x,t) = vect[mesh_matrix.nb_rows()-x-1];
        }
    }


    std::vector<double> solver::solve_tridiag(const matrix& mat, std::vector<double>& d){

        if(mat.nb_cols()!=mat.nb_rows()){
            throw std::domain_error("Not square matrix" );
        }

        if(mat.nb_cols()!=d.size()){
            throw std::domain_error("b of wrong size when solving Ax=b" );
        }
        
        std::size_t n = mat.nb_cols();
        std::vector<double> a(n-1);
        std::vector<double> b(n);
        std::vector<double> c(n-1);

        for(std::size_t i=0; i<n; i++){
            for( std::size_t j=0; j<n;j++){
                if(j==i+1){
                    c[i]=mat(i,j);
                }
                else if(j==i-1){
                    a[j]=mat(i,j);
                }
                else if (j==i){
                    b[i]=mat(i,j);
                }
            }
        }

        std::vector<double> x(n);

        for(int i=1; i<=n-1; i++){
            
            double w = double(a[i-1])/double(b[i-1]);
            b[i] -= w*c[i-1];
            d[i] -= w*d[i-1];
        }
        x[n-1] = d[n-1]/b[n-1];
        for(int i=n-1; i >= 1; i--)
            x[i-1] = (d[i-1]- c[i-1]*x[i])/b[i-1];

		return x;
        
    }





    void solver::compute_price()
    {

        if (m_vol.get_vol().size() == 1 && m_rate.get_rates().size()==1) 
        {
            double a,b,c,d;
            this->init_coeff(a,b,c,d);
            this->transform_coeff(a,b,c,d);

            const int dim = m_msh.get_ndx()-2;

            matrix m_trans_1(dim,dim);

            matrix m_trans_2(dim,dim);

            this->init_matrice_1(m_trans_1,dim,a,b,c,d); 
            this-> init_matrice_2(m_trans_2,dim,a,b,c,d);

            const int ndx = m_msh.get_ndx();
            const int ndt = m_msh.get_ndt();

            matrix mesh_matrix(ndx,ndt);
        
            std::vector<double> final_poff = m_poff(m_msh.get_xaxis());
            std::vector<double> final_vect(ndx-2);

            this->fill_matrix(mesh_matrix,mesh_matrix.nb_cols()-1,final_poff); //Compute complete final cond

            std::copy(final_poff.begin()+1,final_poff.end()-1,final_vect.begin()); // extract only the core part in final_vect

            //std::cout << std::endl << "Mesh à l'initialisation : " << std::endl << mesh_matrix << std::endl;

            std::vector<double> vect(ndx-2); // ATTENTION AU COEFF D
            vect[0] = -m_bd.get_lower_b() * a + d;
            vect[ndx-3] = - m_bd.get_upper_b()* a + d;

            //std::cout << std::endl << "Matrice Trans  1 : " << std::endl << m_trans_1 << std::endl;

            for (int i=(m_msh.get_ndt()-2); i>=0;i--){

                final_vect = m_trans_2 * final_vect;
                final_vect = final_vect + vect;
                final_vect = this->solve_tridiag(m_trans_1,final_vect);
                //std::cout << std::endl << "Payoff hors BC en colonne : " << i+1 << std::endl << final_vect << std::endl;

                std::copy(final_vect.begin(),final_vect.end(),final_poff.begin()+1);

                this->fill_matrix(mesh_matrix,i,final_poff);
            };

            //std::cout << mesh_matrix << std::endl;
            std::cout << "Prix : " << mesh_matrix((ndx-1)/2,0) << std::endl << std::endl;
        
           

        

            


        


        }

        
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
                                dauphine::boundary& bd) // need to define it


{
    Eigen::Matrix<double, n_dx + 1, n_dt + 1> A;
    Eigen::Matrix<double, n_dx + 1, n_dt + 1> B;
    Eigen::Matrix<double, n_dx + 1, n_dt + 1> C;
    Eigen::Matrix<double, n_dx + 1, n_dt + 1> D;

    const double dt = T / n_dt;
    const double dx = (s_max - s_min) / n_dx;

    const double nu1 = dt / dx^2;
    const double nu2 = dt / dx;

    for ( int i=0 ; i < n_dx + 1 ; i++) {
	    for ( int j=0 ; j < n_dt + 1 ; j++ {
		    A(i,j) =  (a(i,j)*nu1 - (1/2) * b(i,j) * nu2) ;
		    B(i,j) =  (c(i,j) * dt - 2 * a(i,j) * nu1) ;
		    C(i,j) =  (a(i,j) * nu1  + (1/2) * b(i,j) * nu2) ;
		    D(i,j) =  d(i,j) * dt ;
	    }
    }


};*/

}

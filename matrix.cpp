#include "matrix.hpp"
#include <math.h>
#include <string>

namespace dauphine {

    std::ostream& operator<<(std::ostream& out, const matrix& m)
    {
        for(std::size_t i = 0; i < m.nb_rows(); ++i)
        {
            for(std::size_t j = 0; j < m.nb_cols(); ++j)
            {
                out << m(i, j) << " ";
            }
            out << std::endl;
        }
        return out;
    }

    matrix operator+(const matrix& lhs, const matrix& rhs)
    {
        matrix tmp(lhs);
        tmp += rhs;
        return tmp;
    }

    matrix operator+(const matrix& rhs, double lhs)
    {
        return rhs + lhs;
    }

    matrix operator+(double lhs, const matrix& rhs)
    {
        return rhs + lhs;
    }

    std::vector<double> operator*(const matrix& lhs, const std::vector<double>& rhs)
    {
        std::vector<double> tmp(lhs.nb_rows());

        if(lhs.nb_cols()!=rhs.size()){
            throw std::domain_error("Wrong vector dimension for multiplication" );
        }
        else{
            for(std::size_t i=0; i<lhs.nb_rows(); i++){
                double temp = 0 ;
                for(std::size_t  j=0;j<lhs.nb_cols(); j++){
                    temp = temp + lhs(i,j)*rhs[j];
                }
                tmp[i]=temp;
            }
        }
        return tmp;
    }

    matrix::matrix(std::size_t nb_rows, std::size_t nb_cols)
        : m_nb_rows(nb_rows),
          m_nb_cols(nb_cols),
          m_data(nb_rows * nb_cols)
    {
    }

    matrix::matrix(const matrix& rhs)
        : m_nb_rows(rhs.m_nb_rows),
          m_nb_cols(rhs.m_nb_cols),
          m_data(rhs.m_data)
    {
    }

    std::size_t matrix::nb_rows() const
    {
        return m_nb_rows;
    }

    std::size_t matrix::nb_cols() const
    {
        return m_nb_cols;
    }

    void matrix::resize(std::size_t nb_rows, std::size_t nb_cols)
    {
        m_nb_rows = nb_rows;
        m_nb_cols = nb_cols;
        m_data.resize(m_nb_rows*m_nb_cols);
    }

    double& matrix::operator()(std::size_t i, std::size_t j)
    {
        return m_data[i * m_nb_cols + j];
    }

    const double& matrix::operator()(std::size_t i, std::size_t j) const
    {
        return m_data[i * m_nb_cols + j];
    }

    matrix& matrix::operator+=(const matrix& rhs)
    {
        std::transform(m_data.begin(), m_data.end(), rhs.m_data.begin(),
                       m_data.begin(), std::plus<double>());
        return *this;
    }

    matrix& matrix::operator-=(const matrix& rhs)
    {
        std::transform(m_data.begin(), m_data.end(), rhs.m_data.begin(),
                       m_data.begin(), std::minus<double>());
        return *this;
    }

    matrix& matrix::operator*=(const matrix& rhs)
    {
        std::transform(m_data.begin(), m_data.end(), rhs.m_data.begin(),
                       m_data.begin(), std::multiplies<double>());
        return *this;
    }

    matrix& matrix::operator/=(const matrix& rhs)
    {
        std::transform(m_data.begin(), m_data.end(), rhs.m_data.begin(),
                       m_data.begin(), std::divides<double>());
        return *this;
    }

    matrix& matrix::operator+=(double rhs)
    {
        std::transform(m_data.begin(), m_data.end(), m_data.begin(),
                       [rhs](double arg) { return arg + rhs; });
        return *this;
    }

    matrix& matrix::operator-=(double rhs)
    {
        std::transform(m_data.begin(), m_data.end(), m_data.begin(),
                       [rhs](double arg) { return arg - rhs; });
        return *this;
    }

    matrix& matrix::operator*=(double rhs)
    {
        std::transform(m_data.begin(), m_data.end(), m_data.begin(),
                       [rhs](double arg) { return arg * rhs; });
        return *this;
    }
    matrix& matrix::operator/=(double rhs)
    {
        std::transform(m_data.begin(), m_data.end(), m_data.begin(),
                       [rhs](double arg) { return arg / rhs; });
        return *this;
    }

}
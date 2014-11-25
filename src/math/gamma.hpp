//**************************************************************************************//
//     Copyright (C) 2014 Malik Kirchner "malik.kirchner@gmx.net"                       //
//                                                                                      //
//     This program is free software: you can redistribute it and/or modify             //
//     it under the terms of the GNU General Public License as published by             //
//     the Free Software Foundation, either version 3 of the License, or                //
//     (at your option) any later version.                                              //
//                                                                                      //
//     This program is distributed in the hope that it will be useful,                  //
//     but WITHOUT ANY WARRANTY; without even the implied warranty of                   //
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                    //
//     GNU General Public License for more details.                                     //
//                                                                                      //
//     You should have received a copy of the GNU General Public License                //
//     along with this program.  If not, see <http://www.gnu.org/licenses/>.            //
//                                                                                      //
//     Dieses Programm ist Freie Software: Sie können es unter den Bedingungen          //
//     der GNU General Public License, wie von der Free Software Foundation,            //
//     Version 3 der Lizenz oder (nach Ihrer Option) jeder späteren                     //
//     veröffentlichten Version, weiterverbreiten und/oder modifizieren.                //
//                                                                                      //
//     Dieses Programm wird in der Hoffnung, dass es nützlich sein wird, aber           //
//     OHNE JEDE GEWÄHRLEISTUNG, bereitgestellt; sogar ohne die implizite               //
//     Gewährleistung der MARKTFÄHIGKEIT oder EIGNUNG FÜR EINEN BESTIMMTEN ZWECK.       //
//     Siehe die GNU General Public License für weitere Details.                        //
//                                                                                      //
//     Sie sollten eine Kopie der GNU General Public License zusammen mit diesem        //
//     Programm erhalten haben. Wenn nicht, siehe <http://www.gnu.org/licenses/>.       //
//                                                                                      //
//**************************************************************************************//


#pragma once

#include <complex>
#include <iostream>
#include <Eigen/Core>

#include <math/pauli.hpp>

namespace math {

constexpr size_t __spinor_dim( const size_t D ) {
    size_t res = D;
    if ( D&1 ) res++;
    return (1 << res/2);
}


template< typename BT, size_t D >
struct __gamma {
    typedef std::integral_constant<size_t, D&1?D-1:D >                   dim;
    typedef std::integral_constant<bool, static_cast<bool>(D&1)>         is_odd;

    typedef std::complex<BT>                                             body_type;
    typedef std::integral_constant<size_t, __spinor_dim(dim::value)>     spinor_dim;
    typedef std::integral_constant<size_t, __spinor_dim(dim::value-2)>   prev_spinor_dim;
    typedef Eigen::Matrix<body_type, spinor_dim::value     , spinor_dim::value     >    gamma_type;
    typedef Eigen::Matrix<body_type, prev_spinor_dim::value, prev_spinor_dim::value>    prev_gamma_type;

    static std::vector<gamma_type> compile() noexcept {
        typedef typename PauliMatrixGenerator<BT>::pauli_type pauli_type;
        const pauli_type sigma1 = PauliMatrixGenerator<BT>::sigma1();
        const pauli_type sigma2 = PauliMatrixGenerator<BT>::sigma2();
        const pauli_type sigma3 = PauliMatrixGenerator<BT>::sigma3();
        const body_type  I{0,1};

        const std::vector<prev_gamma_type> prev_gamma = __gamma<BT, dim::value-2>::compile();
        std::vector<gamma_type> gamma( dim::value + 1 );

        if ( is_odd::value ) {
            gamma[dim::value] *= body_type{0,1};
        } else {
            for ( size_t d = 0; d < dim::value-2; d++) {
                for ( size_t i = 0; i < prev_spinor_dim::value; ++i )
                for ( size_t k = 0; k < prev_spinor_dim::value; ++k ) {
                    gamma[d].template block<2,2>(2*i,2*k) = prev_gamma[d](i,k)*sigma1;
                }
            }

            gamma[dim::value-2] = gamma_type::Zero();
            for ( size_t m = 0; m < spinor_dim::value; m += 2 ) {
                gamma[dim::value-2].template block<2,2>(m,m) = sigma2;
            }

            gamma[dim::value-1] = gamma_type::Zero();
            for ( size_t m = 0; m < spinor_dim::value; m += 2 ) {
                gamma[dim::value-1].template block<2,2>(m,m) = sigma3;
            }

            for ( size_t i = 0; i < prev_spinor_dim::value; i++ )
            for ( size_t k = 0; k < prev_spinor_dim::value; k++ ) {
                gamma[dim::value].template block<2,2>(2*i,2*k) = -prev_gamma[dim::value-2](i,k)*sigma1;
            }

        }

        gamma.shrink_to_fit();

        for ( gamma_type& g : gamma ) {
            for ( size_t m = 0; m < spinor_dim::value; m++ )
                for ( size_t n = 0; n < spinor_dim::value; n++ )
                    if ( fabs(g(m,n)) < 1e-1 ) g(m,n) = body_type{+0.};
        }

        return gamma;
    }
};

template< typename BT >
struct __gamma<BT, 3> {
    typedef std::complex<BT>                                                body_type;
    typedef std::integral_constant<size_t, __spinor_dim(2)>                 spinor_dim;
    typedef Eigen::Matrix<body_type, spinor_dim::value, spinor_dim::value>  gamma_type;

    static std::vector<gamma_type> compile() noexcept {
        typedef typename PauliMatrixGenerator<BT>::pauli_type pauli_type;
        const pauli_type sigma1 = PauliMatrixGenerator<BT>::sigma1();
        const pauli_type sigma2 = PauliMatrixGenerator<BT>::sigma2();
        const pauli_type sigma3 = PauliMatrixGenerator<BT>::sigma3();
        const body_type  I{0,1};

        std::vector<gamma_type> gamma(3+1);

        gamma[0] =   sigma1;
        gamma[1] =   sigma2;
        gamma[2] = I*sigma3;
        gamma.shrink_to_fit();

        for ( gamma_type& g : gamma ) {
            for ( size_t m = 0; m < g.rows(); m++ )
                for ( size_t n = 0; n < g.cols(); n++ )
                    if ( fabs(g(m,n)) < 1e-1 ) g(m,n) = body_type{+0.};
        }

        return gamma;
    }
};

template< typename BT >
struct __gamma<BT, 2> {
    typedef std::complex<BT>                                                body_type;
    typedef std::integral_constant<size_t, __spinor_dim(2)>                 spinor_dim;
    typedef Eigen::Matrix<body_type, spinor_dim::value, spinor_dim::value>  gamma_type;

    static std::vector<gamma_type> compile() noexcept {
        typedef typename PauliMatrixGenerator<BT>::pauli_type pauli_type;
        const pauli_type sigma1 = PauliMatrixGenerator<BT>::sigma1();
        const pauli_type sigma2 = PauliMatrixGenerator<BT>::sigma2();
        const pauli_type sigma3 = PauliMatrixGenerator<BT>::sigma3();
        const body_type  I{0,1};

        std::vector<gamma_type> gamma(2+1);

        gamma[0] =   sigma1;
        gamma[1] =   sigma2;
        gamma[2] = I*sigma3;
        gamma.shrink_to_fit();

        for ( gamma_type& g : gamma ) {
            for ( size_t m = 0; m < g.rows(); m++ )
                for ( size_t n = 0; n < g.cols(); n++ )
                    if ( fabs(g(m,n)) < 1e-1 ) g(m,n) = body_type{+0.};
        }

        return gamma;
    }
};

template< typename BT, size_t D >
class GammaMatrixGenerator {
public:
    typedef std::complex<BT>   body_type;
    typedef BT                 scalar_type;

    typedef std::integral_constant<size_t, D&1?D-1:D >                      dim;
    typedef std::integral_constant<bool, static_cast<bool>(D&1)>            is_odd;

    typedef std::integral_constant<size_t, __spinor_dim(dim::value)>        spinor_dim;
    typedef typename PauliMatrixGenerator<BT>::pauli_type                   pauli_type;
    typedef Eigen::Matrix<body_type,spinor_dim::value,spinor_dim::value>    gamma_type;

private:
    const std::vector<gamma_type> _gamma;

public:

    GammaMatrixGenerator() :
        _gamma( __gamma<BT, D>::compile() )
    {

    }

    gamma_type operator[]( size_t k ) const {
        assert( k < D );
        return _gamma[k];
    }

    gamma_type chiral() const {
        assert( !is_odd::value/*, "There is no chirality operator in odd dimensions, "
                                       "since the square charge conjugation operator is purely imaginary." */);
        return _gamma[D];
    }


    void print() {
        for ( gamma_type g : _gamma ) {
            std::cout << g << std::endl << std::endl;
        }
    }

};



}

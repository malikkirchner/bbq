/*****************************************************************************************
==========================================================================================
                               ____      ____      ___
                            U | __")u U | __")u   / " \
                             \|  _ \/  \|  _ \/  | |"| |
                              | |_) |   | |_) | /| |_| |\
                              |____/    |____/  U \__\_\u
                             _|| \\_   _|| \\_     \\//
                            (__) (__) (__) (__)   (_(__)


                               -- The Fermion grill --
==========================================================================================
*****************************************************************************************/



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
#include <Eigen/Sparse>

#include <math/pauli.hpp>

namespace math {

//! matrix dimension of the spinor rep. in D space time dimensions
constexpr size_t __spinor_dim( const size_t D ) {
    return D&1 ? (1 << (D-1)/2) : (1 << D/2);
}

template< typename BT, size_t D >
class __gamma {
private:
    //! canonic space time dimension of the spinor representation
    typedef std::integral_constant<size_t, D&1?D-1:D >                                  canonic_dim;
    //! Is the actual space time dimension odd?
    typedef std::integral_constant<bool  , static_cast<bool>(D&1)>                      is_odd;
    //! complex number type (body)
    typedef std::complex<BT>                                                            body_type;
    //! scalar type, complex numbers are based on
    typedef std::complex<BT>                                                            scalar_type;
    //! matrix dimension of the spinor rep. in D space time dimensions
    typedef std::integral_constant<size_t, __spinor_dim(canonic_dim::value)>            spinor_dim;
    //! matrix dimension of the spinor rep. in D-2 space time dimensions
    typedef std::integral_constant<size_t,
                    __spinor_dim(canonic_dim::value-(is_odd::value?0:2))>               prev_spinor_dim;
    //! type of the gamma matrices in D dimensions
    typedef Eigen::Matrix<body_type, spinor_dim::value     , spinor_dim::value     >    gamma_type;
    //! type of the gamma matrices in D-2 dimensions
    typedef Eigen::Matrix<body_type, prev_spinor_dim::value, prev_spinor_dim::value>    prev_gamma_type;

    static void fill( const std::true_type t, std::vector<gamma_type>& gamma,
                      const std::vector<prev_gamma_type>& prev_gamma )
    {
        const body_type  I{0,1};
        gamma.assign( prev_gamma.begin(), prev_gamma.end() );
        gamma[canonic_dim::value] *= body_type{0,1};
    }

    static void fill( const std::false_type t, std::vector<gamma_type>& gamma,
                      const std::vector<prev_gamma_type>& prev_gamma )
    {
        const body_type  I{0,1};
        typedef typename PauliMatrixGenerator<BT>::pauli_type pauli_type;
        const pauli_type sigma1 = PauliMatrixGenerator<BT>::sigma1();
        const pauli_type sigma2 = PauliMatrixGenerator<BT>::sigma2();
        const pauli_type sigma3 = PauliMatrixGenerator<BT>::sigma3();

        const size_t dim      = __spinor_dim(canonic_dim::value);
        const size_t prev_dim = __spinor_dim(canonic_dim::value-2);

        for ( size_t d = 0; d < canonic_dim::value-2; ++d ) {
            for ( size_t i = 0; i < prev_dim; ++i )
            for ( size_t k = 0; k < prev_dim; ++k ) {
                gamma[d].template block<2,2>(2*i,2*k) = prev_gamma[d](i,k)*sigma1;
            }
        }

        gamma[canonic_dim::value-2] = gamma_type::Zero();
        for ( size_t m = 0; m < dim; m += 2 ) {
            gamma[canonic_dim::value-2].template block<2,2>(m,m) = sigma2;
        }

        gamma[canonic_dim::value-1] = gamma_type::Zero();
        for ( size_t m = 0; m < dim; m += 2 ) {
            gamma[canonic_dim::value-1].template block<2,2>(m,m) = sigma3;
        }

        for ( size_t i = 0; i < prev_dim; ++i )
        for ( size_t k = 0; k < prev_dim; ++k ) {
            gamma[canonic_dim::value].template block<2,2>(2*i,2*k) = -prev_gamma[canonic_dim::value-2](i,k)*sigma1;
        }
    }

public:
    static std::vector<gamma_type> compile() noexcept {
        const std::vector<prev_gamma_type> prev_gamma = __gamma<BT, canonic_dim::value-(is_odd::value?0:2)>::compile();
        std::vector<gamma_type> gamma( canonic_dim::value + 1 );

        fill( typename is_odd::type(), gamma, prev_gamma );

        gamma.shrink_to_fit();

        for ( gamma_type& g : gamma ) {
            for ( size_t m = 0; m < spinor_dim::value; ++m )
                for ( size_t n = 0; n < spinor_dim::value; ++n )
                    if ( fabs(g(m,n)) < 1e-1 ) g(m,n) = body_type{+0.};
        }

        return gamma;
    }
};

template< typename BT >
class __gamma<BT, 3> {
private:
    //! canonic space time dimension of the spinor representation
    typedef std::integral_constant<size_t, 3&1?3-1:3 >                          canonic_dim;
    //! Is the actual space time dimension odd?
    typedef std::integral_constant<bool  , static_cast<bool>(3&1)>              is_odd;
    //! complex number type (body)
    typedef std::complex<BT>                                                    body_type;
    //! scalar type, complex numbers are based on
    typedef std::complex<BT>                                                    scalar_type;
    //! matrix dimension of the spinor rep. in D space time dimensions
    typedef std::integral_constant<size_t, __spinor_dim(canonic_dim::value)>    spinor_dim;
    //! type of the gamma matrices in D dimensions
    typedef Eigen::Matrix<body_type, spinor_dim::value, spinor_dim::value >     gamma_type;

public:
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
            for ( size_t m = 0; m < spinor_dim::value; ++m )
                for ( size_t n = 0; n < spinor_dim::value; ++n )
                    if ( fabs(g(m,n)) < 1e-1 ) g(m,n) = body_type{+0.};
        }

        return gamma;
    }
};

template< typename BT >
class __gamma<BT, 2> {
private:
    //! canonic space time dimension of the spinor representation
    typedef std::integral_constant<size_t, 2&1?2-1:2 >                          canonic_dim;
    //! Is the actual space time dimension odd?
    typedef std::integral_constant<bool  , static_cast<bool>(2&1)>              is_odd;
    //! complex number type (body)
    typedef std::complex<BT>                                                    body_type;
    //! scalar type, complex numbers are based on
    typedef std::complex<BT>                                                    scalar_type;
    //! matrix dimension of the spinor rep. in D space time dimensions
    typedef std::integral_constant<size_t, __spinor_dim(canonic_dim::value)>    spinor_dim;
    //! type of the gamma matrices in D dimensions
    typedef Eigen::Matrix<body_type, spinor_dim::value, spinor_dim::value >     gamma_type;

public:
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
            for ( size_t m = 0; m < spinor_dim::value; ++m )
                for ( size_t n = 0; n < spinor_dim::value; ++n )
                    if ( fabs(g(m,n)) < 1e-1 ) g(m,n) = body_type{+0.};
        }

        return gamma;
    }
};


/*!**************************************************************************************
 * @class  GammaMatrixGenerator
 * @author Malik Kirchner <malik.kirchner@gmx.net>
 *
 * @brief  Generates Dirac gamma matrices in arbitrary dimension.
 *
 * http://en.wikipedia.org/wiki/Higher-dimensional_gamma_matrices
 ****************************************************************************************/
template< typename BT, size_t D >
class GammaMatrixGenerator {
public:
    //! canonic space time dimension of the spinor representation
    typedef std::integral_constant<size_t, D&1?D-1:D >                                  canonic_dim;
    //! Is the actual space time dimension odd?
    typedef std::integral_constant<bool  , static_cast<bool>(D&1)>                      is_odd;
    //! complex number type (body)
    typedef std::complex<BT>                                                            body_type;
    //! scalar type, complex numbers are based on
    typedef std::complex<BT>                                                            scalar_type;
    //! matrix dimension of the spinor rep. in D space time dimensions
    typedef std::integral_constant<size_t, __spinor_dim(canonic_dim::value)>            spinor_dim;
    //! type of the gamma matrices in D dimensions
    typedef Eigen::Matrix<body_type, spinor_dim::value, spinor_dim::value >             gamma_type;
    //! type of the gamma matrices in D-2 dimensions in sparse rep.
    typedef Eigen::SparseMatrix<body_type>                                              sparse_gamma_type;

private:
    const std::vector<gamma_type>        _gamma;
    const std::vector<sparse_gamma_type> _gamma_sparse;

    std::vector<sparse_gamma_type> makeSparse( std::vector<gamma_type> in ) {
        std::vector<sparse_gamma_type> out( in.size() );

        for ( size_t k = 0; k < in.size(); ++k )
            out[k] = in[k].sparseView(1e-2);

        return out;
    }

public:

    GammaMatrixGenerator() noexcept :
        _gamma( __gamma<BT, D>::compile() ),
        _gamma_sparse( makeSparse(_gamma) )
    {

    }

    const gamma_type& operator[]( const size_t k ) const noexcept {
        assert( k < D );
        return _gamma[k];
    }

    const gamma_type& chiral() const noexcept {
        assert( !is_odd::value );
        return _gamma[D];
    }

    void print() {
        for ( gamma_type g : _gamma ) {
            std::cout << g << std::endl << std::endl;
        }
    }

};



}

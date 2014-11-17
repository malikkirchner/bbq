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


#include <math/sun.hpp>

namespace math {



/*!**************************************************************************************
 *
 * @class  suBase
 * @author Malik Kirchner <malik.kirchner@gmx.net>
 * @brief  Base in su(N)-algebra.
 *
 * suBase is a template for a class representing a su(N)-algebra base. It
 * provides generation of an orthogonal base in su(N), especially the Pauli
 * matrices for su(2) and the Gell-Mann matrices for su(3). Use suBase to
 * project arbitrary su(N) elements or create a su(N) element as linear
 * combination of the orthogonal base.
 * http://en.wikipedia.org/wiki/Generalizations_of_Pauli_matrices
 *
 ****************************************************************************************/
template< typename BT, size_t N >
class suNBase {
protected:

    constexpr size_t addr( const size_t m, const size_t n ) const noexcept {
        return N*m+n;
    }

public:
    typedef typename su<BT,N>::body_type body_type;

    BT          traces  [N*N];
    su<BT,N>    base    [N*N];
    struct Coefficients { BT c[N*N]; };

    suNBase() {
        // construct a base in a su(N) (eg. Gell-Mann Matrices for su(3))
        for (size_t k = 0; k < N; k++) {
            for (size_t j = 0; j < N; j++) {
                if ( k < j ) {
                    base[addr(k,j)].setZero();
                    base[addr(k,j)](k,j) = body_type(1.0,0.0);
                    base[addr(k,j)](j,k) = body_type(1.0,0.0);
                }

                if ( k == j ) {
                    if ( k == 0 ) {
                        base[addr(k,j)].setZero();
                        for ( size_t zz = 0; zz < N; zz++)
                            base[addr(k,k)](zz,zz) = body_type(1.0,0.0);
                    } else if ( k > 0 ) {
                        base[addr(k,j)].setZero();
                        for ( size_t zz = 0; zz < k; zz++)
                            base[addr(k,k)](zz,zz) = sqrt(2.0/(BT)(k*(k+1)))*body_type(1.0,0.0);
                        base[addr(k,k)](k,k)       = sqrt(2.0/(BT)(k*(k+1)))*body_type(-(BT)k,0.0);
                    }
                }

                if ( k > j ) {
                    base[addr(k,j)].setZero();
                    base[addr(k,j)](k,j) = body_type(0.0,+1.0);
                    base[addr(k,j)](j,k) = body_type(0.0,-1.0);
                }
            }
        }

        // pre-calc traces of squared base matrices for normalization in project() ...
        for (size_t k = 0; k<N*N; k++) {
            traces[k] = base[k].norm2();
        }
    }


    const Coefficients project(const su< BT, N >& rhs) const {
        su< BT, N >     buf;
        Coefficients    coeff;
        BT              aux     = 1.0;

        for (size_t k = 0; k < N*N; k++) {
            // multiply
            buf = rhs*base[k];
            // trace
            aux = 2.0/traces[k];
            coeff.c[k] = 0.0;
            for (size_t j = 0; j < N; j++)
                coeff.c[k] += aux*buf(j,j).imag();
        }

        return coeff;
    }

    const su< BT, N > span(const Coefficients& rhs) const {
        su< BT, N > buf;
        buf.setZero();
        // span
        const body_type I = body_type(0.0,5.0);
        for (size_t k = 0; k < N*N; k++) {
            buf += I*rhs.c[k]*base[k];
        }
        return buf;
    }

};



}

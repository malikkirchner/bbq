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


#include <vector>
#include <random>
#include <complex>
#include <math/matrixfactory.hpp>
#include <math/su.hpp>
#include <math/subase.hpp>

namespace math {



/*!**************************************************************************************
 * @class  SUFactory
 * @author Malik Kirchner <malik.kirchner@gmx.net>
 *
 * @brief  SU(N) matrix generator.
 ****************************************************************************************/
template<typename BT, size_t N>
class SUFactory : public MatrixFactory< SU<BT, N> > {
protected:

    std::mt19937 uni;

    inline const BT inverseDist(const BT x) {
        BT y                   = (BT)0.0;
        BT a;
        BT a0;
        const BT c0 = 1.8171205928321397;
        const BT c1 = 0.4;
        const BT c2 = 0.22641786849563159;
        const BT c3 = 0.16613673991608136;
        const BT c4 = 0.13785528756957333;
        const BT c5 = 0.12291798578358654;
        const BT c6 = 0.11490033055470908;
        const BT c7 = 0.11105453691182769;
        const BT c8 = 0.11004173524517248;
        a  = pow(x, 1.0/3.0);   //a^1/3
        a0 = a*a;
        y += c0*a;
        a *= a0;                //a^3/3
        y += c1*a;
        a *= a0;                //a^5/3
        y += c2*a;
        a *= a0;                //a^7/3
        y += c3*a;
        a *= a0;                //a^9/3
        y += c4*a;
        a *= a0;                //a^11/3
        y += c5*a;
        a *= a0;                //a^13/3
        y += c6*a;
        a *= a0;                //a^15/3
        y += c7*a;
        a *= a0;                //a^17/3
        y += c8*a;

        return y;
    }

public:
    typedef SU<BT, N>   matrix_type;

    virtual matrix_type generateRandom() noexcept final {
        BT rangemod = 0.5;
        const std::complex<BT> I = std::complex<BT>(0,1);

        matrix_type res;

        BT cos_alpha     = (BT)0.0;
        BT cos_phi       = (BT)0.0;
        BT cos_theta     = (BT)0.0;

        BT sin_alpha     = (BT)0.0;
        BT sin_phi       = (BT)0.0;
        BT sin_theta     = (BT)0.0;

        BT phi           = (BT)0.0;
        BT alpha         = (BT)0.0;

        if ( N == 1 ) {
            cos_alpha = static_cast< BT >(uni());
            res(0,0).real() = cos_alpha;
            res(0,0).imag() = sqrt(1.0-cos_alpha*cos_alpha);
        } else {
            math::suBase<BT,2>                         base;
            math::su<BT,2>                             G;
            math::SU<BT,N>                             A;
            typename math::suBase<BT,2>::Coefficients  coeff;

            res = identity();
            for ( unsigned i = 0; i < N-1; i++ ) {
                for ( unsigned j = i+1; j < N; j++ ) {
                    alpha     = inverseDist((BT)0.3930385334019695*static_cast< BT >(uni()));
                    phi       = (BT)2.0* M_PI*static_cast< BT >(uni()-0.5);
                    cos_theta = (BT)2.0*static_cast< BT >(uni())-1.0;

                    cos_alpha = (BT)1.0-(BT)(1.0+(2.0*(BT)(rand()%2)-1.0)*cos(alpha))*rangemod;
                    sin_alpha = (BT)sqrt(1.0-cos_alpha*cos_alpha);
                    sin_theta = (BT)sqrt(1.0-cos_theta*cos_theta);
                    cos_phi   = (BT)cos(phi);
                    sin_phi   = (BT)sin(phi);

                    coeff.c[0] = (BT)sin_alpha*sin_theta*cos_phi;
                    coeff.c[1] = (BT)sin_alpha*sin_theta*sin_phi;
                    coeff.c[2] = (BT)sin_alpha*cos_theta;

                    G  = cos_alpha;	              ; G -= I*coeff.c[0]*base.base[1];
                    G -= I*coeff.c[1]*base.base[2]; G -= I*coeff.c[2]*base.base[3];

                    A = res;
                    for ( unsigned k = 0; k < N; k++ ) {
                        res(i,k)  = G(0,0) * A(i,k);
                        res(i,k) += G(0,1) * A(j,k);
                        res(j,k)  = G(1,0) * A(i,k);
                        res(j,k) += G(1,1) * A(j,k);
                    }
                }
            }
        }

        return res;
    }

    virtual matrix_type identity() noexcept final {
        return matrix_type::Identity();
    }

    virtual matrix_type zero() noexcept final {
        matrix_type res;
        res.setZero();
        return res;
    }

};

}

/*****************************************************************************************
==========================================================================================
                               ____      ____      ___
                            U | __")u U | __")u   / " \
                             \|  _ \/  \|  _ \/  | |"| |
                              | |_) |   | |_) | /| |_| |\
                              |____/    |____/  U \__\_\u
                             _|| \\_   _|| \\_     \\//
                            (__) (__) (__) (__)   (_(__)


                               -- The Fermion Grill --
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

#include <Eigen/Core>
#include <complex>
#include <iostream>

namespace math {

/*!**************************************************************************************
 * @class  PauliMatrixGenerator
 * @author Malik Kirchner <malik.kirchner@gmx.net>
 *
 * @brief  Generates Pauli matrices.
 ****************************************************************************************/
template <typename BT> class PauliMatrixGenerator {
public:
    using body_type   = std::complex<BT>;
    using scalar_type = BT;
    using pauli_type  = Eigen::Matrix<body_type, 2, 2>;

    //! \return Pauli matrix \f$\sigma_1\f$
    static pauli_type sigma1() noexcept {
        pauli_type m;
        m(0, 0) = body_type{0, 0};
        m(0, 1) = body_type{1, 0};
        m(1, 0) = body_type{1, 0};
        m(1, 1) = body_type{0, 0};
        return m;
    }

    //! \return Pauli matrix \f$\sigma_2\f$
    static pauli_type sigma2() noexcept {
        pauli_type m;
        m(0, 0) = body_type{0, 0};
        m(0, 1) = body_type{0, -1};
        m(1, 0) = body_type{0, 1};
        m(1, 1) = body_type{0, 0};
        return m;
    }

    //! \return Pauli matrix \f$\sigma_1\f$
    static pauli_type sigma3() noexcept {
        pauli_type m;
        m(0, 0) = body_type{1, 0};
        m(0, 1) = body_type{0, 0};
        m(1, 0) = body_type{0, 0};
        m(1, 1) = body_type{-1, 0};
        return m;
    }
};
}

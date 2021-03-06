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

#include <operators/wilsondirac.hpp>

namespace operators {

/*!**************************************************************************************
 * @class  OverlapDirac
 * @author Malik Kirchner <malik.kirchner@gmx.net>
 *
 * @brief  Overlap-Dirac operator.
 ****************************************************************************************/
template <typename fermion_traits, typename gauge_traits> class OverlapDirac : public Operator {
public:
    using body_type     = typename Operator<fermion_traits, gauge_traits>::body_type;
    using scalar_type   = typename Operator<fermion_traits, gauge_traits>::scalar_type;
    using fermion_field = typename Operator<fermion_traits, gauge_traits>::fermion_field;
    using gauge_field   = typename Operator<fermion_traits, gauge_traits>::gauge_field;
    using gauge_field   = typename Operator<fermion_traits, gauge_traits>::gauge_field;
    using base_operator = WilsonDirac<fermion_traits, gauge_traits>;

private:
    base_operator _baseOperator;

public:
    virtual fermion_field apply(const fermion_field& phi, const gauge_field& U) const final {
        fermion_field res(phi);

        return res;
    }
};
}

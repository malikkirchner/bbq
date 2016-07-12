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

#include <field/basefield.hpp>
#include <math/sufactory.hpp>

namespace field {

template <class MatrixType, class LatticeType, field_periodicity Periodicity>
struct gauge_field_traits : public field_traits<MatrixType, LatticeType, true, Periodicity> {};

template <class MatrixType, class LatticeType>
struct periodic_gauge_field_traits
    : public gauge_field_traits<MatrixType, LatticeType, FP_PERIODIC> {};

template <class MatrixType, class LatticeType>
struct anti_periodic_gauge_field_traits
    : public gauge_field_traits<MatrixType, LatticeType, FP_ANTI_PERIODIC> {};

template <class MatrixType, class LatticeType>
struct finite_gauge_field_traits : public gauge_field_traits<MatrixType, LatticeType, FP_FINITE> {};

/*!**************************************************************************************
 * @class  GaugeField
 * @author Malik Kirchner <malik.kirchner@gmx.net>
 *
 * @brief  Gauge field over lattice links.
 ****************************************************************************************/
template <class Traits> class GaugeField : public BaseField<Traits> {
public:
    using base_type    = BaseField<Traits>;
    using traits       = Traits;
    using gauge_type   = typename traits::matrix_type;
    using lattice_type = typename traits::lattice_type;

private:
protected:
public:
    GaugeField(const lattice_type& lattice)
        : base_type(lattice) {}

    GaugeField(const GaugeField& other)
        : base_type(other) {}

    inline gauge_type& operator()(const std::size_t x, const std::size_t mu) noexcept {
        return base_type::operator[](x * 4 + mu);
    }

    inline gauge_type const& operator()(const std::size_t x, const std::size_t mu) const noexcept {
        return base_type::operator[](x * 4 + mu);
    }

    inline std::size_t numLinks() const noexcept { return base_type::_volume * base_type::_dim; }
};
}

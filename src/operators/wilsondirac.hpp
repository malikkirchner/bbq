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

#include <type_traits>

#include <field/neighbours.hpp>
#include <math/gamma.hpp>
#include <operators/operator.hpp>

namespace operators {

/*!**************************************************************************************
 * @class  WilsonDirac
 * @author Malik Kirchner <malik.kirchner@gmx.net>
 *
 * @brief  Wilson-Dirac operator.
 ****************************************************************************************/
template <typename fermion_traits, typename gauge_traits>
class WilsonDirac : public Operator<fermion_traits, gauge_traits> {
public:
    using base_type     = Operator<fermion_traits, gauge_traits>;
    using body_type     = typename base_type::body_type;
    using scalar_type   = typename base_type::scalar_type;
    using fermion_field = typename base_type::fermion_field;
    using gauge_field   = typename base_type::gauge_field;
    using lattice_type  = typename base_type::lattice_type;
    using lattice_dim   = typename fermion_traits::lattice_dim;

private:
    const lattice_type _lattice;

    scalar_type _mass;
    scalar_type _r;

    using periodicity = typename fermion_traits::periodicity;
    const field::Neighbours<lattice_type, periodicity::value> neighbours;

    math::GammaMatrixGenerator<scalar_type, lattice_dim::value> gamma;

public:
    WilsonDirac(lattice_type lattice_)
        : _lattice(lattice_)
        , _mass(1.)
        , _r(1.)
        , neighbours(lattice_) {}

    WilsonDirac(lattice_type lattice_, scalar_type m_, scalar_type r_)
        : _lattice(lattice_)
        , _mass(m_)
        , _r(r_)
        , neighbours(lattice_) {}

    virtual fermion_field apply(const fermion_field& phi, const gauge_field& U) const final {
        assert(phi.lattice() == _lattice);
        assert(phi.lattice() == U.lattice());

        const std::size_t volume = _lattice.volume();
        const std::size_t dim    = _lattice.dim();
        fermion_field     res(_lattice);

        for (std::size_t k = 0; k < volume; k++) {
            res[k] = phi[k];
            res[k] *= (4.0 * _r + _mass);
        }

        for (std::size_t m = 0; m < volume; m++) {
            for (std::size_t n = 0; n < volume; n++) {
                for (std::size_t mu = 0; mu < dim; mu++) {}
            }
        }

        return res;
    }

    scalar_type mass() const noexcept { return _mass; }
    scalar_type r() const noexcept { return _r; }

    void mass(const scalar_type other) noexcept { _mass = other; }
    void r(const scalar_type other) noexcept { _r = other; }

    lattice_type lattice() const noexcept { return _lattice; }
};
}

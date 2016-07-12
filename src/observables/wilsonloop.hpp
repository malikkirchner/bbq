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

#include <field/neighbours.hpp>
#include <observables/observable.hpp>

namespace observable {

/*!**************************************************************************************
 * @class  Observable
 * @author Malik Kirchner <malik.kirchner@gmx.net>
 *
 * @brief  Wilson loop observable
 *
 * http://en.wikipedia.org/wiki/Wilson_loop
 ****************************************************************************************/
template <typename Gauge, std::size_t W = 1u, std::size_t H = 1u>
class WilsonLoop : public Observable {
public:
    using gauge_field     = Gauge;
    using gauge_type      = typename gauge_field::gauge_type;
    using lattice_type    = typename gauge_field::traits::lattice_type;
    using body_type       = typename gauge_field::traits::body_type;
    using periodicity     = typename gauge_field::traits::periodicity;
    using neighbours_type = field::Neighbours<lattice_type, periodicity::value>;

    using width  = std::integral_constant<std::size_t, W>;
    using height = std::integral_constant<std::size_t, H>;

private:
    static constexpr auto loop(const gauge_field& gauge, const neighbours_type& neighbours,
                               const std::size_t x, const std::size_t mu, const std::size_t nu,
                               std::true_type is_plaquette) noexcept {
        return gauge(x, nu).adjoint() * gauge(neighbours(x, nu).fwd, mu).adjoint() *
               gauge(neighbours(x, mu).fwd, nu) * gauge(x, mu);
    }

    static constexpr gauge_type loop(const gauge_field& gauge, const neighbours_type& neighbours,
                                     const std::size_t x, const std::size_t mu,
                                     const std::size_t nu, std::false_type is_plaquette) noexcept {
        gauge_type L;
        L.setIdentity();

        std::size_t xx = x;

        for (std::size_t w = 0; w < width::value; w++) {
            L *= gauge(xx, mu);
            xx = neighbours(xx, mu).fwd;
        }

        for (std::size_t h = 0; h < height::value; h++) {
            L *= gauge(xx, nu);
            xx = neighbours(xx, nu).fwd;
        }

        for (std::size_t w = 0; w < width::value; w++) {
            L *= gauge(xx, mu).adjoint();
            xx = neighbours(xx, mu).bwd;
        }

        for (std::size_t h = 0; h < height::value; h++) {
            L *= gauge(xx, nu).adjoint();
            xx = neighbours(xx, nu).bwd;
        }

        return L;
    }

public:
    static body_type eval(const gauge_field& gauge, const neighbours_type& neighbours) noexcept {
        const lattice_type lattice = gauge.lattice();

        const std::size_t volume = lattice.volume();
        const std::size_t dim    = lattice_type::lattice_dim::value;

        gauge_type wl;
        wl.setZero();

        using is_plaquette = typename std::integral_constant<bool, (width::value == 1u) &&
                                                                       (height::value == 1u)>::type;

        for (std::size_t x = 0; x < volume; ++x) {
            for (std::size_t mu = 0; mu < dim; ++mu)
                for (std::size_t nu = 0; nu < mu; ++nu) {
                    wl += loop(gauge, neighbours, x, mu, nu, is_plaquette());
                }
        }

        return wl.trace();
    }
};
}

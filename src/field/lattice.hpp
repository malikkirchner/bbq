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

#include <cassert>
#include <math/mod.hpp>
#include <type_traits>

namespace field {

/*!**************************************************************************************
 * @class  Lattice
 * @author Malik Kirchner <malik.kirchner@gmx.net>
 *
 * @brief  Space time lattice in arbitrary dimension.
 ****************************************************************************************/
template <std::size_t D, typename fast_mod = std::true_type> class Lattice {
public:
    template <typename T = long> struct Index {
        using index_type = T;

        index_type c[D];
        inline index_type& operator[](const index_type k) noexcept { return c[k]; }
        inline index_type const& operator[](const index_type k) const noexcept { return c[k]; }

        Index()             = default;
        Index(const Index&) = default;
        Index(Index&&)      = default;
        Index& operator=(const Index&) = default;
        Index& operator=(Index&&) = default;
    };

    using lattice_dim = std::integral_constant<std::size_t, D>;

private:
    template <typename T> constexpr std::size_t __volume__(const Index<T>& idx) const noexcept {
        std::size_t res = idx[0];
        for (std::size_t k = 1; k < D; k++) res *= idx[k];
        return res;
    }

protected:
    const Index<std::size_t> _dimension;
    const std::size_t        _volume;

public:
    Lattice(const Index<std::size_t>& dim) noexcept : _dimension(dim), _volume(__volume__(dim)) {
        static_assert(D > 1, "Lattice dimension must be at least two.");
        assert(_volume > 0);
    }

    Lattice(const Lattice&) = default;
    Lattice(Lattice&&)      = default;
    Lattice& operator=(const Lattice&) = default;
    Lattice& operator=(Lattice&&) = default;

    constexpr bool operator==(const Lattice& other) const noexcept {
        return (other._volume == _volume);
    }

    constexpr std::size_t addr(const Index<std::size_t>& idx) const noexcept {
        std::size_t res = idx[0];

        for (std::size_t k = 1; k < D; k++) {
            res *= _dimension[k];
            res += idx[k];
        }

        return res;
    }

    inline std::size_t addr_mod(const Index<long>& idx) const noexcept {
        std::size_t res = 0;

        if (fast_mod::value) {

            res = math::mod(idx[0], static_cast<long>(_dimension[0]));

            for (std::size_t k = 1; k < D; k++) {
                res *= _dimension[k];
                res += math::mod(idx[k], static_cast<long>(_dimension[k]));
            }

        } else {

            res = ((idx[0] % _dimension[0]) + _dimension[0]) % _dimension[0];

            for (std::size_t k = 1; k < D; k++) {
                res *= _dimension[k];
                res += ((idx[k] % _dimension[k]) + _dimension[k]) % _dimension[k];
            }
        }

        return res;
    }

    constexpr std::size_t dim() const noexcept { return lattice_dim::value; }
    constexpr std::size_t volume() const noexcept { return _volume; }
    constexpr std::size_t dimension(const std::size_t k) const noexcept { return _dimension[k]; }
};
}

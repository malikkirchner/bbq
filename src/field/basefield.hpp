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
#include <cstring>
#include <type_traits>
#include <util/memory.hpp>

namespace field {

enum field_periodicity { FP_PERIODIC = 1, FP_ANTI_PERIODIC = 2, FP_FINITE = 4 };

template <class MatrixType, class LatticeType, bool IsLink, field_periodicity Periodicity>
struct field_traits {
    using matrix_type  = MatrixType;
    using lattice_type = LatticeType;
    using body_type    = typename MatrixType::body_type;
    using scalar_type  = typename MatrixType::scalar_type;
    using is_link      = std::integral_constant<bool, IsLink>;
    using periodicity  = std::integral_constant<field_periodicity, Periodicity>;
    using lattice_dim  = typename LatticeType::lattice_dim;
};

/*!**************************************************************************************
 * @class  GaugeField
 * @author Malik Kirchner <malik.kirchner@gmx.net>
 *
 * @brief  Basic field over sites or links of a lattice.
 ****************************************************************************************/
template <class Traits> class BaseField {
public:
    using traits       = Traits;
    using matrix_type  = typename traits::matrix_type;
    using lattice_type = typename traits::lattice_type;

protected:
    matrix_type*       _data;
    const lattice_type _lattice;
    const std::size_t  _dim;
    const std::size_t  _volume;

    void allocate() {
        const std::size_t len = traits::is_link::value ? _dim * _volume : _volume;
        _data                 = new matrix_type[len];
    }

    void deallocate() { safe_array_delete(_data); }

    void cloneData(typename traits::matrix_type* other) {
        const std::size_t len = traits::is_link::value ? _dim * _volume : _volume;
        for (std::size_t k = 0; k < len; k++) _data[k] = other[k];
    }

public:
    BaseField(lattice_type lattice)
        : _data(nullptr)
        , _lattice(lattice)
        , _dim(_lattice.dim())
        , _volume(_lattice.volume()) {
        allocate();
    }

    BaseField(const BaseField& other)
        : _data(nullptr)
        , _lattice(other._lattice)
        , _dim(_lattice.dim())
        , _volume(_lattice.volume()) {
        allocate();
        cloneData(other._data);
    }

    BaseField& operator=(const BaseField& other) {
        assert(_lattice == other._lattice);

        if (!_data) allocate();
        cloneData(other._data);

        return *this;
    }

    inline matrix_type& operator[](const std::size_t k) noexcept { return _data[k]; }

    inline matrix_type const& operator[](const std::size_t k) const noexcept { return _data[k]; }

    virtual ~BaseField() { deallocate(); }

    constexpr long                dim() const noexcept { return _dim; }
    constexpr long                volume() const noexcept { return _volume; }
    constexpr lattice_type const& lattice() const noexcept { return _lattice; }
};
}

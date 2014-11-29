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
#include <cstring>
#include <cassert>
#include <util/memory.hpp>


namespace field {
    
enum field_periodicity {
    FP_PERIODIC = 1,
    FP_ANTI_PERIODIC = 2,
    FP_FINITE = 4
};
    
template<class MatrixType, class LatticeType, bool IsLink, field_periodicity Periodicity>
struct field_traits {
    typedef MatrixType                          matrix_type;
    typedef LatticeType                         lattice_type;
    typedef typename MatrixType::body_type      body_type;
    typedef typename MatrixType::scalar_type    scalar_type;
    typedef std::integral_constant< bool, IsLink >                    is_link;
    typedef std::integral_constant< field_periodicity, Periodicity >  periodicity;
    typedef typename LatticeType::lattice_dim                         lattice_dim;
};
    
        
/*!**************************************************************************************
 * @class  GaugeField
 * @author Malik Kirchner <malik.kirchner@gmx.net>
 *
 * @brief  Basic field over sites or links of a lattice.
 ****************************************************************************************/
template< class Traits >
class BaseField {
public:
    typedef Traits                          traits;
    typedef typename traits::matrix_type    matrix_type;
    typedef typename traits::lattice_type   lattice_type;
    
protected:
    matrix_type*        _data;
    const lattice_type  _lattice;
    const size_t        _dim;
    const size_t        _volume;
    
    void allocate() {
        const size_t len = traits::is_link::value ? _dim*_volume : _volume;
        _data = new matrix_type [ len ];
    }

    void deallocate() {
        safe_array_delete( _data );
    }

    void cloneData( typename traits::matrix_type *other ) {
        const size_t len = traits::is_link::value ? _dim*_volume : _volume;
        for ( size_t k = 0; k < len; k++ ) _data[k] = other[k];
    }

public:    
    
    BaseField( lattice_type lattice ) :
        _data( NULL ), _lattice( lattice ),
        _dim( _lattice.dim() ), _volume( _lattice.volume() )
    {
        allocate();
    }
    
    BaseField( const BaseField& other ) :
        _data( NULL ), _lattice( other._lattice ),
        _dim( _lattice.dim() ), _volume( _lattice.volume() )
    {
        allocate();
        cloneData( other._data );
    }

    BaseField& operator = ( const BaseField& other ) {
        assert( _lattice == other._lattice );

        if ( !_data ) allocate();
        cloneData( other._data );

        return *this;
    }

    constexpr matrix_type& operator[] ( const size_t k ) noexcept {
        return _data[k];
    }

    constexpr matrix_type const & operator[] ( const size_t k ) const noexcept {
        return _data[k];
    }


    virtual ~BaseField() {
        deallocate();
    }
    
    constexpr long dim()    const noexcept { return _dim; }
    constexpr long volume() const noexcept { return _volume; }
    constexpr lattice_type const & lattice() const noexcept { return _lattice; }
};
    
}

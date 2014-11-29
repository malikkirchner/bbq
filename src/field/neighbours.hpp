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

#include <cstring>
#include <util/memory.hpp>
#include <field/basefield.hpp>

namespace field {


/*!**************************************************************************************
 * @class  Neighbours
 * @author Malik Kirchner <malik.kirchner@gmx.net>
 *
 * @brief  Generates next neighbour field.
 ****************************************************************************************/
template< class LT, field_periodicity Periodicity >
class Neighbours {
public:

    struct Entry {
        size_t fwd;
        size_t bwd;
    };

    typedef LT                                          lattice_type;
    typedef typename lattice_type::template Index<long> index_type;

protected:

    const lattice_type  _lattice;
    Entry*              _data;
    const size_t        _dim;
    const size_t        _volume;

    inline size_t addr( const index_type& idx, const long d ) const noexcept {
        return _lattice.addr_mod(idx)*_dim + d;
    }

    inline size_t addr( const size_t idx, const size_t d ) const noexcept {
        return idx*_dim + d;
    }

    void compile() noexcept {
        index_type idx;
        size_t     mu = 0;

        compile( idx, mu );
    }

    void compile( index_type& idx, const size_t mu ) noexcept {
        const long L = _lattice.dimension(mu);

        if ( mu < _dim - 1 ) {
            for ( idx[mu] = 0; idx[mu] < L; ++idx[mu] ) {
                compile( idx, mu + 1 );
            }
        } else {
            for ( idx[mu] = 0; idx[mu] < L; ++idx[mu] ) {
                index_type nb_idx = idx;
                for ( size_t d = 0; d < _dim; ++d ) {
                    ++nb_idx[d];
                    _data[ addr( idx, d ) ].fwd = _lattice.addr_mod( nb_idx );
                    nb_idx[d] -= 2;
                    _data[ addr( idx, d ) ].bwd = _lattice.addr_mod( nb_idx );
                }
            }
        }
    }

public:

    Neighbours( lattice_type lattice_ ) :
        _lattice( lattice_ ), _data(NULL),
        _dim(_lattice.dim()), _volume(_lattice.volume())
    {
        _data = new Entry [ _dim*_volume ];
        compile();
    }

    Neighbours( const Neighbours& other ) :
        _lattice( other._lattice ), _data(NULL),
        _dim(_lattice.dim()), _volume(_lattice.volume())
    {
        _data = new Entry [ _dim*_volume ];
        memcpy( _data, other._data, _dim*_volume*sizeof(Entry) );
    }

    Neighbours operator = ( const Neighbours& other ) {
        assert( _lattice == other._lattice );
        _data = new Entry [ _dim*_volume ];
        memcpy( _data, other._data, _dim*_volume*sizeof(Entry) );
    }

    Neighbours( Neighbours&& other ) = default;
    Neighbours& operator = ( Neighbours&& other ) = default;

    virtual ~Neighbours() {
        safe_array_delete( _data );
    }

    const Entry& operator () ( const size_t     k, const size_t d ) const noexcept { return _data[ addr(k,d) ]; }
    const Entry& operator () ( const index_type k, const size_t d ) const noexcept { return _data[ addr(k,d) ]; }

    lattice_type lattice() const noexcept { return _lattice; }
};

}

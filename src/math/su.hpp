﻿//**************************************************************************************//
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

#include <complex>
#include <Eigen/Core>

namespace math {


/*!**************************************************************************************
 * @author Malik Kirchner <malik.kirchner@gmx.net>
 *
 ****************************************************************************************/
template< typename BT, size_t N >
class SU {
private:
    Eigen::Matrix< std::complex<BT>, N, N > data;

public:
    typedef std::complex<BT>                        body_type;
    typedef Eigen::Matrix< std::complex<BT>, N, N > matrix_type;

    constexpr SU<BT, N> dagger() {
        Eigen::MatrixBase< matrix_type >  mb(data);
        mb.adjointInPlace();
        return mb;
    }

    void setZero() { data.setZero(); }

    constexpr body_type& operator()( const size_t m, const size_t n ) noexcept { return data(m,n); }
    constexpr body_type const & operator()( const size_t m, const size_t n ) const noexcept { return data(m,n); }
};


/*!**************************************************************************************
 * @author Malik Kirchner <malik.kirchner@gmx.net>
 *
 ****************************************************************************************/
template< typename BT, size_t N >
class su  {
private:
    Eigen::Matrix< std::complex<BT>, N, N > data;

public:
    typedef std::complex<BT>                        body_type;
    typedef Eigen::Matrix< std::complex<BT>, N, N > matrix_type;

    constexpr su<BT, N> dagger() {
        Eigen::MatrixBase< matrix_type > mb(data);
        mb.adjointInPlace();
        return su{ data=mb };
    }

    void setZero() { data.setZero(); }

    constexpr body_type& operator()( const size_t m, const size_t n ) noexcept { return data(m,n); }
    constexpr body_type const & operator()( const size_t m, const size_t n ) const noexcept { return data(m,n); }

    constexpr BT norm2 () {
        const Eigen::MatrixBase< matrix_type > mb(data);
        return (mb*data).trace().real();
    }

    std::ostream& print( std::ostream& os ) {
        os << data;
        return os;
    }

};



}
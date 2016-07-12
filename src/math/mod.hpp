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

namespace math {

inline int mod(const int idx, const int d) noexcept {
    const int dm1 = d - 1;
    int       aux = idx;
    aux &= dm1;
    aux += d;
    aux &= dm1;
    return aux;
}

inline long mod(const long idx, const long d) noexcept {
    const long dm1 = d - 1;
    long       aux = idx;
    aux &= dm1;
    aux += d;
    aux &= dm1;
    return aux;
}

inline unsigned mod(const unsigned idx, const unsigned d) noexcept {
    unsigned aux = idx;
    aux &= d - 1;
    return aux;
}

inline unsigned long mod(const unsigned long idx, const unsigned long d) noexcept {
    unsigned long aux = idx;
    aux &= d - 1;
    return aux;
}
}

/*
 * BENCHMARK 29.11.2014
 *
 * Lattice 32x32x32x32
 *
 * Intel(R) Core(TM) i7 CPU 920 @ 2.67GHz
 *
 * CLang 3.5.0
 * GCC   4.9.2
 *
 * CXX_FLAGS=-march=native -fomit-frame-pointer -pipe -Wall -pedantic -g0 -O2
 * -fpic -std=c++14 -DNDEBUG
 *
 *              |      GCC      |     CLang
 * -------------|---------------|---------------
 * fast mod     |     23.2s     |      2.4s
 * stl  mod     |     91.4s     |     91.5s
 *
 */

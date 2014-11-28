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


#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE "lattice-test"

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <chrono>
#include <field/lattice.hpp>



#define LAT_ADDR(a,b,c,d) (((( (a%32+32)%32 )*32+( (b%32+32)%32 ))*32+( (c%32+32)%32 ))*32+( (d%32+32)%32 ))

BOOST_AUTO_TEST_CASE(LatticePerformanceTest) {
    using namespace field;

    Lattice<4> lat{{{32,32,32,32}}};

    double  stl_mod;
    double  fast_mod;
    long x0 = 0;
    long x1 = 0;

    {
        auto start = std::chrono::high_resolution_clock::now();

        for ( long a = -32; a < 500*32; a++ )
        for ( long b = -32; b < 32; b++ )
        for ( long c = -32; c < 32; c++ )
        for ( long d = -32; d < 32; d++ ) {
            x0 += lat.addr_mod({{a,b,c,d}});
        }
        auto stop   = std::chrono::high_resolution_clock::now();
        fast_mod    = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start).count();

        std::cout << "fast_mod: " << fast_mod << std::endl;
    }

    {
        auto start = std::chrono::high_resolution_clock::now();

        for ( long a = -32; a < 500*32; a++ )
        for ( long b = -32; b < 32; b++ )
        for ( long c = -32; c < 32; c++ )
        for ( long d = -32; d < 32; d++ ) {
            x1 += LAT_ADDR(a,b,c,d);
        }
        auto stop   = std::chrono::high_resolution_clock::now();
        stl_mod     = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start).count();

        std::cout << "stl_mod : " << stl_mod << std::endl;
    }

    std::cout << "stl_mod/fast_mod : " << (double)stl_mod/(double)fast_mod << std::endl;

    BOOST_CHECK_EQUAL( x0, x1 );
}

BOOST_AUTO_TEST_CASE(LatticeTest) {
    using namespace field;

    Lattice<4> lat{{{2,8,32,128}}};

    const long a = lat.addr({{1,3,7,127}});
    const long b = ((1*8+3)*32+7)*128+127;
    const long c = lat.addr_mod({{1,3+8,7-32,127}});

    BOOST_CHECK_EQUAL( a, b );
    BOOST_CHECK_EQUAL( b, c );
}


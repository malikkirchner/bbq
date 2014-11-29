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


#include <iostream>
#include <field/lattice.hpp>
#include <field/gaugefield.hpp>
#include <field/fermionfield.hpp>
#include <math/spinor.hpp>
#include <math/su.hpp>
#include <math/subase.hpp>
#include <operators/wilsondirac.hpp>
#include <observables/wilsonloop.hpp>

/*!**************************************************************************************
 * @author Malik Kirchner <malik.kirchner@gmx.net>
 *
 ****************************************************************************************/
int main ( int argc, char** argv ) {
    int EXIT_CODE = 0;

    field::Lattice<4> lattice{{{8,8,8,8}}};

    typedef field::Lattice<4>           lattice_type;
    typedef math::Spinor<double, 4, 3>  spinor_type;
    typedef math::SU<double, 3>         gauge_type;
    typedef field::periodic_fermion_field_traits< spinor_type, lattice_type >   fermion_traits;
    typedef field::periodic_gauge_field_traits< gauge_type, lattice_type >      gauge_traits;
    typedef field::GaugeField< gauge_traits >       gauge_field;
    typedef field::FermionField< fermion_traits >   fermion_field;

    fermion_field   phi(lattice);
    gauge_field     U(lattice);
    operators::WilsonDirac< fermion_traits, gauge_traits > wilson(lattice);

    phi = wilson.apply( phi, U );


    field::Neighbours< gauge_traits::lattice_type, gauge_traits::periodicity::value > neighbours( lattice );
    observable::WilsonLoop<gauge_field>::eval( U, neighbours );

    using namespace math;

    math::suBase<double,3> sub;
    for ( size_t k = 0; k < 3*3; k++ )
        sub.base[k].print(std::cout) << std::endl << std::endl;

    math::GammaMatrixGenerator<double, 4> gamma;

//    gamma.print();

    return EXIT_CODE;
}

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


#include <field/basefield.hpp>

namespace field {


template< class MatrixType, class LatticeType, field_periodicity Periodicity >
struct fermion_field_traits : public field_traits< MatrixType, LatticeType, false, Periodicity >{};

template< class MatrixType, class LatticeType >
struct periodic_fermion_field_traits : public fermion_field_traits< MatrixType, LatticeType, FP_PERIODIC >{};

template< class MatrixType, class LatticeType >
struct anti_periodic_fermion_field_traits : public fermion_field_traits< MatrixType, LatticeType, FP_ANTI_PERIODIC >{};

template< class MatrixType, class LatticeType >
struct finite_fermion_field_traits : public fermion_field_traits< MatrixType, LatticeType, FP_FINITE >{};
    
/*!**************************************************************************************    
 * @author Malik Kirchner <malik.kirchner@gmx.net>
 * 
 ****************************************************************************************/    
template< class Traits >
class FermionField : public BaseField< Traits > {
public:
    typedef BaseField< Traits >          base_type;
    typedef Traits                       traits;
    typedef typename traits::matrix_type spinor_type;

private:
    
protected:
        
public:
    
    FermionField( const typename traits::lattice_type& lattice ) : base_type( lattice ) 
    {
        
    }
    
};
    
}
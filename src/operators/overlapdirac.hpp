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

#include <operators/wilsondirac.hpp>

namespace operators {


template< typename fermion_traits, typename gauge_traits >
class OverlapDirac : public Operator {
public:
    typedef typename Operator< fermion_traits, gauge_traits >::body_type     body_type;
    typedef typename Operator< fermion_traits, gauge_traits >::scalar_type   scalar_type;
    typedef typename Operator< fermion_traits, gauge_traits >::fermion_field fermion_field;
    typedef typename Operator< fermion_traits, gauge_traits >::gauge_field   gauge_field;
    typedef typename Operator< fermion_traits, gauge_traits >::gauge_field   gauge_field;
    typedef WilsonDirac< fermion_traits, gauge_traits >                      base_operator;

private:
    base_operator _baseOperator;

public:

    virtual fermion_field apply( const fermion_field& phi, const gauge_field& U ) const final {
        fermion_field res( phi );



        return res;
    }

};

}

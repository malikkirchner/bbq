##**************************************************************************************##
##     Copyright (C) 2014 Malik Kirchner "malik.kirchner@gmx.net"                       ##
##                                                                                      ##
##     This program is free software: you can redistribute it and/or modify             ##
##     it under the terms of the GNU General Public License as published by             ##
##     the Free Software Foundation, either version 3 of the License, or                ##
##     (at your option) any later version.                                              ##
##                                                                                      ##
##     This program is distributed in the hope that it will be useful,                  ##
##     but WITHOUT ANY WARRANTY; without even the implied warranty of                   ##
##     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                    ##
##     GNU General Public License for more details.                                     ##
##                                                                                      ##
##     You should have received a copy of the GNU General Public License                ##
##     along with this program.  If not, see <http:##www.gnu.org/licenses/>.            ##
##                                                                                      ##
##     Dieses Programm ist Freie Software: Sie können es unter den Bedingungen          ##
##     der GNU General Public License, wie von der Free Software Foundation,            ##
##     Version 3 der Lizenz oder (nach Ihrer Option) jeder späteren                     ##
##     veröffentlichten Version, weiterverbreiten und/oder modifizieren.                ##
##                                                                                      ##
##     Dieses Programm wird in der Hoffnung, dass es nützlich sein wird, aber           ##
##     OHNE JEDE GEWÄHRLEISTUNG, bereitgestellt; sogar ohne die implizite               ##
##     Gewährleistung der MARKTFÄHIGKEIT oder EIGNUNG FÜR EINEN BESTIMMTEN ZWECK.       ##
##     Siehe die GNU General Public License für weitere Details.                        ##
##                                                                                      ##
##     Sie sollten eine Kopie der GNU General Public License zusammen mit diesem        ##
##     Programm erhalten haben. Wenn nicht, siehe <http:##www.gnu.org/licenses/>.       ##
##                                                                                      ##
##**************************************************************************************##

add_definitions(-DBOOST_TEST_DYN_LINK) 
# clang fix
add_definitions( "-D__extern_always_inline=inline" )


add_executable( latticetest latticetest.cpp )
target_link_libraries( latticetest ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} )
set_target_properties( latticetest PROPERTIES RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_BINARY_DIR}/test-bin)

enable_testing()
add_test( NAME latticetest WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test-bin COMMAND latticetest )

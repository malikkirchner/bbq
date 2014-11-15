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


cmake_minimum_required(VERSION 2.8)

set(CMAKE_C_COMPILER   "/usr/bin/clang"  )
set(CMAKE_CXX_COMPILER "/usr/bin/clang++")


project(test)

set(CMAKE_CXX_FLAGS "-march=native -pipe -Wall -pedantic -std=c++14 -g1 -O2")

add_executable(bbq src/main.cpp)

install(TARGETS bbq RUNTIME DESTINATION bin)

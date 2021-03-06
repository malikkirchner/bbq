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

#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE "gamma-test"

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <math/gamma.hpp>

template <typename gamma_type>
bool checkEta(size_t mu, size_t nu, gamma_type commutator, double tol = 1e-4) {
    bool res = true;

    const size_t L = static_cast<size_t>(commutator.cols());
    for (size_t i = 0; i < L; ++i)
        for (size_t k = 0; k < L; ++k) {
            if ((i != k) || (mu != nu)) res &= fabs(commutator(i, k)) < tol;
            if ((i == k) && (mu == nu)) res &= 2. - fabs(commutator(i, k)) < tol;
        }

    return res;
}

template <typename BT, size_t D>
bool checkSignature(const math::GammaMatrixGenerator<BT, D>& gamma, BT tol = 1e-4) {
    bool res = true;

    const size_t L = math::GammaMatrixGenerator<BT, D>::spinor_dim::value;
    for (size_t d = 0; d < D; ++d) {
        typename math::GammaMatrixGenerator<BT, D>::gamma_type id = gamma[d] * gamma[d];
        for (size_t i = 0; i < L; ++i)
            for (size_t k = 0; k < L; ++k) {
                if (i != k) res &= fabs(id(i, k)) < tol;
                if (i == k) res &= 1. - fabs(id(i, k)) < tol;
            }
    }

    return res;
}

template <typename BT, size_t D>
bool checkChirality(const math::GammaMatrixGenerator<BT, D>& gamma, BT tol = 1e-4) {
    if (math::GammaMatrixGenerator<BT, D>::is_odd::value) return true;

    bool res = true;

    typedef typename math::GammaMatrixGenerator<BT, D>::gamma_type gamma_type;
    gamma_type ch = pow(std::complex<BT>{0., 1.}, D / 2 - 1) * gamma[0];

    for (size_t d = 1; d < D; ++d) ch *= gamma[d];

    const size_t L = math::GammaMatrixGenerator<BT, D>::spinor_dim::value;
    for (size_t i = 0; i < L; ++i)
        for (size_t k = 0; k < L; ++k) { res &= fabs(ch(i, k) - gamma.chiral()(i, k)) < tol; }

    const gamma_type unit = 2. * gamma.chiral() * gamma.chiral();
    BOOST_CHECK_MESSAGE(checkEta(D, D, unit), "Chirality matrix is not self inverse!");

    for (size_t d = 0; d < D; ++d) {
        const gamma_type com = gamma[d] * ch + ch * gamma[d];

        bool aux = checkEta(D, d, com);
        BOOST_CHECK_MESSAGE(aux, "Chirality matrix does not anti-commute with gamma matrices!");
        if (!aux) std::cout << com << std::endl;
    }

    std::cout << "\n--------------------------------------------------\n";
    std::cout << "expected chirality matrix\n";
    std::cout << ch << std::endl;
    std::cout << "actual chirality matrix\n";
    std::cout << gamma.chiral() << std::endl;

    return res;
}

template <typename BT, size_t D> void testGamma() {
    using namespace math;

    GammaMatrixGenerator<BT, D>                              gamma;
    typedef typename GammaMatrixGenerator<BT, D>::gamma_type gamma_type;

    std::cout << "\n";
    for (size_t mu = 0; mu < D; mu++)
        for (size_t nu = 0; nu < D; nu++) {
            const gamma_type com = gamma[mu] * gamma[nu] + gamma[nu] * gamma[mu];

            std::cout << "\n--------------------------------------------------\n";
            std::cout << "(mu,nu) = (" << mu << "," << nu << ")\n";
            std::cout << com << std::endl;

            BOOST_CHECK_MESSAGE(checkEta(mu, nu, com),
                                "Euclidean anti-commutation relation is broken.");
        }

    BOOST_CHECK_MESSAGE(checkChirality(gamma), "Chirality matrix is invalid.");

    BOOST_CHECK_MESSAGE(checkSignature(gamma), "Spacetime signature is invalid.");
}

BOOST_AUTO_TEST_CASE(GammaTest) {
    testGamma<double, 2>();
    testGamma<double, 3>();
    testGamma<double, 4>();
    testGamma<double, 5>();
    testGamma<double, 6>();
}

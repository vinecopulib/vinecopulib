//
// Created by Vatter Thibault on 14/12/16.
//


#include <common/include/joe_bicop.h>
#include <iostream>

int main(int __attribute__((unused)) argc, __attribute__((unused)) char *argv[]) {

    JoeBicop a;
    std::cout << "A joe copula...\n";
    std::cout << "Family : " << a.get_family() << "\n";
    std::cout << "Parameter : " << a.get_parameters() << "\n";
    std::cout << "Rotation : " << a.get_rotation() << "\n";
    std::cout << "Number of parameters : " << a.calculate_npars() << "\n \n";

    MatXd U = MatXd::Zero(2,2);
    U(0,0) = 0.3;
    U(0,1) = 0.4;
    U(1,0) = 0.8;
    U(1,1) = 0.6;
    std::cout << "U = \n" << U << "\n";
    std::cout << "pdf = \n" << a.pdf(U) << "\n";
    std::cout << "hfunc1 = \n" << a.hfunc1(U) << "\n";
    std::cout << "hfunc2 = \n" << a.hfunc2(U) << "\n \n";

    double theta = 3;
    JoeBicop b(theta);
    std::cout << "Another joe copula...\n";
    std::cout << "Family : " << b.get_family() << "\n";
    std::cout << "Parameter : " << b.get_parameters() << "\n";
    std::cout << "Rotation : " << b.get_rotation() << "\n";
    std::cout << "Number of parameters : " << b.calculate_npars() << "\n \n";

    std::cout << "U = \n" << U << "\n";
    std::cout << "generator for u1 = \n" << b.generator(U.col(0)) << "\n";
    std::cout << "generator_inv for u1 = \n" << b.generator_inv(U.col(0)) << "\n";
    std::cout << "generator_derivative for u1 = \n" << b.generator_derivative(U.col(0)) << "\n";
    std::cout << "generator_derivative2 for u1 = \n" << b.generator_derivative2(U.col(0)) << "\n";
    std::cout << "pdf = \n" << b.pdf(U) << "\n";
    std::cout << "hfunc1 = \n" << b.hfunc1(U) << "\n";
    std::cout << "hfunc2 = \n" << b.hfunc2(U) << "\n";
    std::cout << "hinv1 = \n" << b.hinv1(U) << "\n";
    std::cout << "hinv2 = \n" << b.hinv2(U) << "\n \n";

    return 0;
}

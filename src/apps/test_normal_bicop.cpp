//
// Created by Vatter Thibault on 14/12/16.
//


#include <common/include/normal_bicop.h>
#include <iostream>

int main(int __unused argc, char __unused *argv[]) {

    NormalBicop a;
    std::cout << "A normal copula...\n";
    std::cout << "Family : " << a.get_family() << "\n";
    std::cout << "Parameter : " << a.get_parameters() << "\n";
    std::cout << "Number of parameters : " << a.calculate_npars() << "\n \n";

    /*MatXd U = MatXd::Zero(6,2);
    U(0,0) = 0.3;
    U(0,1) = 1.0;
    U(1,0) = 0.7;
    U(1,1) = 1.0;
    U(2,0) = 1.0;
    U(2,1) = 0.3;
    U(3,0) = 1.0;
    U(3,1) = 0.7;
    U(4,0) = 0.3;
    U(4,1) = 0.7;
    U(5,0) = 0.7;
    U(5,1) = 0.3;*/
    MatXd U = MatXd::Zero(2,2);
    U(0,0) = 0.3;
    U(0,1) = 0.4;
    U(1,0) = 0.8;
    U(1,1) = 0.6;
    std::cout << "U = \n" << U << "\n";
    std::cout << "pdf = \n" << a.pdf(U) << "\n";
    std::cout << "hfunc1 = \n" << a.hfunc1(U) << "\n";
    std::cout << "hfunc2 = \n" << a.hfunc2(U) << "\n \n";

    double rho = 0.5;
    NormalBicop b(rho);
    std::cout << "Another normal copula...\n";
    std::cout << "Family : " << b.get_family() << "\n";
    std::cout << "Parameter : " << b.get_parameters() << "\n";
    std::cout << "Number of parameters : " << b.calculate_npars() << "\n \n";

    std::cout << "U = \n" << U << "\n";
    std::cout << "pdf = \n" << b.pdf(U) << "\n";
    std::cout << "hfunc1 = \n" << b.hfunc1(U) << "\n";
    std::cout << "hfunc2 = \n" << b.hfunc2(U) << "\n";
    std::cout << "hinv1 = \n" << b.hinv1(U) << "\n";
    std::cout << "hinv2 = \n" << b.hinv2(U) << "\n";

    std::cout << "A random sample : " << b.simulate(10) << "\n";
    return 0;
}

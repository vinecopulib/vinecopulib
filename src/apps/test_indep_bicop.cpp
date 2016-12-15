//
// Created by Vatter Thibault on 14/12/16.
//


#include <common/include/indep_bicop.h>
#include <iostream>

int main(int __unused argc, char __unused *argv[]) {

    IndepBicop a;
    std::cout << "An independance copula...\n";
    std::cout << "Family : " << a.get_family() << "\n";
    std::cout << "Parameter : " << a.get_parameters() << "\n";
    std::cout << "Number of parameters : " << a.calculate_npars() << "\n \n";

    MatXd U = MatXd::Zero(2,2);
    U(0,0) = 0.3;
    U(0,1) = 0.7;
    U(1,0) = 0.7;
    U(1,1) = 0.3;
    std::cout << "U = \n" << U << "\n";
    std::cout << "pdf = \n" << a.pdf(U) << "\n";
    std::cout << "hfunc1 = \n" << a.hfunc1(U) << "\n";
    std::cout << "hfunc2 = \n" << a.hfunc2(U) << "\n";
    return 0;
}

//
// Created by Vatter Thibault on 13/12/16.
//

#include <common/include/parbicop.hpp>
#include <iostream>

int main(int argc, char *argv[]) {

    // Create a dummy object
    ParBiCop a;

    std::cout << "Family : " << a.getFamily() << "\n";
    std::cout << "Parameter : " << a.getPar() << "\n";
    std::cout << "Parameter 2 : " << a.getPar2() << "\n";
    std::cout << "Number of parameters : " << a.getNpars() << "\n";
    return 0;
}

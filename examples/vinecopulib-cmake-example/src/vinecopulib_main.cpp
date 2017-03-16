#include <bicop_class.hpp>
#include <iostream>

using namespace vinecopulib;

int main() {

    // print all available families
	auto family = bicop_families::all;
	std::cout << "Available families : ";
	for (auto fam : family) {
		std::cout << get_family_name(fam) << " ";
	}
	std::cout << "\n";

	// create Clayton copula and simulate from the model
	auto model = Bicop::create(BicopFamily::clayton, 0, Eigen::VectorXd::Constant(1,1));
    auto data = model->simulate(1e3);
    std::cout << "Created Model | family: " << model->get_family_name() << ", ";
	std::cout << "rotation: " <<  model->get_rotation() << ", ";
    std::cout << "parameter: " << model->get_parameters() << std::endl;

	// fit and select family
	auto fitted = Bicop::select(data);
	std::cout << "Created Model | family: " << fitted->get_family_name() << ", ";
	std::cout << "rotation: " <<  fitted->get_rotation() << ", ";
	std::cout << "parameter: " << fitted->get_parameters() << std::endl;

	return 0;
}

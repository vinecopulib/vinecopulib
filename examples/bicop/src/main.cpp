#include <iostream>
#include <vinecopulib.hpp>

using namespace vinecopulib;

int
main()
{

  // print all available families
  auto family = bicop_families::all;
  std::cout << "Available families : ";
  for (auto fam : family) {
    std::cout << get_family_name(fam) << " ";
  }
  std::cout << std::endl;

  // create Gumbel copula and simulate from the model
  auto model = Bicop(BicopFamily::gumbel, 0, Eigen::VectorXd::Constant(1, 2));
  auto data = model.simulate(2e3);
  std::cout << "Created Model | " << model.str() << std::endl;

  // select family and fit
  auto fitted = Bicop(data);
  std::cout << "Fitted Model (sample size = 2e3) | " << fitted.str()
            << std::endl;

  // check the default fit controls
  auto controls = FitControlsBicop();
  std::cout << controls.str() << std::endl;

  return 0;
}

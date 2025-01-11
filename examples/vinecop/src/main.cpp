// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <iostream>
#include <vinecopulib.hpp>

using namespace vinecopulib;

int
main()
{
  // specify the dimension of the model
  int d = 3;

  // instantiate a three dimensional D-vine with independence copulas
  Vinecop default_model(d);

  // alternatively, instantiate a std::vector<std::vector<Bicop>> object
  auto pair_copulas = Vinecop::make_pair_copula_store(d);

  // specify the pair copulas
  auto par = Eigen::VectorXd::Constant(1, 3.0);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 270, par);
    }
  }

  // specify a structure matrix
  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(3, 3);
  mat << 1, 1, 1, 2, 2, 0, 3, 0, 0;

  // instantiate a custom model using pair_copulas and
  Vinecop custom_model(pair_copulas, mat);

  // simulate data
  Eigen::MatrixXd data = custom_model.simulate(1e3);

  // instantiate a D-vine and select the families
  Vinecop fitted(d);
  fitted.select_families(data);

  // alternatively, instantiate a new object by
  // selecting the structure along with the families
  Vinecop fitted2(data);

  // check the default fit controls
  auto controls = FitControlsVinecop();
  std::cout << controls.str() << std::endl;

  return 0;
}

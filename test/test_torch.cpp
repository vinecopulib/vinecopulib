// Copyright © 2016-2020 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <iostream>
#include <torch/torch.h>

int
main()
{
  torch::Tensor tensor = torch::rand({ 2, 3 });
  std::cout << tensor << std::endl;

  std::cout << tensor.sum(1) << std::endl;
  torch::Tensor cuda_tensor = tensor.cuda();
  std::cout << cuda_tensor.sum(1) << std::endl;
}

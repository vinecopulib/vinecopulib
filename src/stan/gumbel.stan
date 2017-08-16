#include <bicop/class.hpp>
functions {
  real gumbel_loglik(matrix[N,2] u, real theta) {
      Bicop model(BicopFamily::gumbel, 0, VecXd::Constant(1, theta));
      return model.loglik(u);
    }

}
data {
  int<lower=0> N;
  matrix[N,2] u;
}
parameters {
  real<lower=1> theta;
}
model {
  target += gumbel_loglik(u, theta);
}

#include <Eigen/Dense>
#include <algorithm>
#include <chrono>

// Computes the median of an Eigen::VectorXd
double
median(Eigen::VectorXd& v)
{
  // Sort the vector in-place
  std::sort(v.data(), v.data() + v.size());
  // Compute and return the median
  return v.size() % 2 == 0 ? v.segment((v.size() - 2) / 2, 2).mean()
                           : v(v.size() / 2);
}

// Overload for const Eigen::VectorXd
double
median(const Eigen::VectorXd& v)
{
  Eigen::VectorXd v_copy = v; // Make a copy of the vector
  return median(v_copy);      // Delegate to the non-const version
}

double
std_dev(const Eigen::VectorXd& v)
{
  double mean = v.mean();
  double sum = (v.array() - mean).square().sum();
  return sqrt(sum / (v.size() - 1));
}

// Benchmarking function with seed support using Eigen::VectorXi
template<typename Function>
Eigen::VectorXd
benchmark_func(Function func, const Eigen::VectorXi& seeds)
{
  unsigned repeats = seeds.size();
  Eigen::VectorXd times = Eigen::VectorXd::Zero(repeats);

  for (unsigned i = 0; i < repeats; ++i) {
    auto start = std::chrono::high_resolution_clock::now();
    func(seeds(i)); // Pass the current seed to the function
    auto end = std::chrono::high_resolution_clock::now();
    times(i) = std::chrono::duration<double, std::milli>(end - start).count();
  }

  return times;
}

Eigen::VectorXd
benchmark_stats(const Eigen::VectorXd& times)
{
  Eigen::VectorXd stats = Eigen::VectorXd::Zero(3);
  stats(0) = times.mean();
  stats(1) = median(times);
  stats(2) = std_dev(times);
  return stats;
}
#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

class MPC final {
public:
  MPC() = delete;

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  static std::vector<double> solve(const Eigen::VectorXd& state, const Eigen::VectorXd& coeffs,
    std::vector<double>& mpcX, std::vector<double>& mpcY);
};

#endif /* MPC_H */
